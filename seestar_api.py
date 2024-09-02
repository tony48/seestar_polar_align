import requests
import json
import time
import tzlocal
from astropy.coordinates import SkyCoord, FK5
from datetime import datetime, timezone
from astropy.time import Time
from astropy import units as u

def send_request(action, parameters):
    ip_address = "localhost"
    port = 5555
    device_num = 1
    base_url = f"http://{ip_address}:{port}/api/v1/telescope/{device_num}/action"

    headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"}
    body = {"Action":action,"Parameters":json.dumps(parameters),"ClientID":1,"ClientTransactionID":999}
    r = requests.put(base_url,body,headers=headers)
    return r

def start_up_sequence():
    send_request("action_start_up_sequence",{"lat":0, "lon":0})

def auto_focus():
    send_request("method_sync",{"method":"start_auto_focuse"})

def stop_auto_focus():
    send_request("method_sync",{"method":"stop_auto_focuse"})

def is_seestar_connected():
    r = send_request("method_sync",{"method":"get_device_state"})
    return r.ok

def move_ra_90(right):
    total_move = 90
    subtotal_move = 0
    if right:
        angle = 0
    else:
        angle = 180
    while subtotal_move < total_move:
        cur_move = min(45, total_move - subtotal_move)
        move_time = cur_move * 20.0 / 90.0
        move_joystick(1000, angle, round(move_time))
        subtotal_move += cur_move
        time.sleep(move_time + 1)

def move_joystick(speed, angle, duration):
    send_request("method_sync", {"method":"scope_speed_move","params":{"speed":speed,"angle":angle,"dur_sec":duration}})

def set_horizontal_calibration(enable):
    send_request("method_sync",{"method":"set_setting", "params":{"auto_3ppa_calib": enable}})

def set_exposure(exposure_stack, exposure_live=2000):
    send_request("method_sync",{"method":"set_setting", "params":{"exp_ms": {"stack_l": exposure_stack,"continuous": exposure_live}}})

def set_lp_filter(state):
    send_request("method_sync",{"method":"set_setting","params":{"stack_lenhance":state}})

def get_wheel_position():
    r = send_request("method_sync",{"method":"get_wheel_position"})
    return r.json()["Value"]["result"]

def get_setting():
    r = send_request("method_sync",{"method":"get_setting"})
    return r.json()

def start_stack(gain=80):
    send_request("start_stack",{"gain":gain, "restart": True})

def stop_stack():
    send_request("method_sync",{"method":"iscope_stop_view","params":{"stage":"Stack"}})

def set_target_name(name):
    send_request("method_sync",{"method":"set_sequence_setting", "params":[{"group_name":name}]})

def get_last_image_url():
    r = send_request("get_last_image",{"is_subframe":True, "is_thumb":False})
    return r.json()["Value"]["url"]

def delete_image(path):
    send_request("method_sync",{"method":"delete_image", "params":{"name":[f"MyWorks/{path}"]}})

def goto(target_name, ra, dec, is_j2000=True):
    # ra and dec are in hms and dms respectively
    parameters = {"target_name":target_name, "ra":ra, "dec":dec, "is_j2000":is_j2000}
    send_request("goto_target",parameters)

def goto_without_platesolve(target_name, ra, dec, is_j2000=True):
    # ra and dec are in hms and dms respectively
    set_target_name(target_name)
    target_coords = parse_coordinate(is_j2000, ra, dec)
    parameters = {"method":"scope_goto","params":[target_coords.ra.hour,target_coords.dec.degree]}
    send_request("method_sync",parameters)

def wait_for_goto():
    print("Waiting for goto to complete")
    while is_goto():
        time.sleep(2)

def wait_for_goto_ra_dec(target_ra, target_dec, is_j2000=True):
    print("Waiting for goto ra and dec")
    target_coords = parse_coordinate(is_j2000, target_ra, target_dec)
    ra = target_coords.ra.hour
    dec = target_coords.dec.degree
    coords = get_equ_coords()
    while (coords["ra"] - ra) ** 2 + (coords["dec"] - dec) ** 2 > 0.01:
        time.sleep(2)
        coords = get_equ_coords()
    time.sleep(1)

def get_equ_coords():
    r = send_request("method_sync", {"method":"scope_get_equ_coord"})
    return r.json()["Value"]["result"]

def is_goto():
    r = send_request("method_sync",{"method": "get_app_state"})
    return r.json()["Value"]["result"]["auto_goto"]["is_working"]

def set_dithering(enable, distance=45, interval=5):
    parameters = {"method":"set_setting", "params":{"stack_dither": {"pix": distance, "interval": interval, "enable": enable}}}
    send_request("method_sync",parameters)

# When changing exposure time, the seestar will shoot dark frames again when start_stack is called
# This function will shoot dark frames if it is necessary
# Call this after set_exposure
def do_darks():
    start_stack()
    time.sleep(2)
    if get_wheel_position() == 0:
        print("Doing darks...")
        while get_wheel_position() == 0:
                time.sleep(2)
        print("Darks are done")
    else:
        print("Darks have already been taken")
    stop_stack()

def parse_coordinate(is_j2000, in_ra, in_dec):
    _fk5 = FK5(equinox=Time(Time(datetime.now(timezone.utc), scale='utc').jd, format="jd", scale="utc"))
    if is_j2000:
        coord_frame = 'icrs'
    else:
        coord_frame = _fk5
    if isinstance(in_ra, str):
        result = SkyCoord(ra=in_ra, dec=in_dec, unit=(u.hourangle, u.deg), frame=coord_frame)
    else:
        result = SkyCoord(ra=in_ra*u.hour, dec=in_dec*u.deg, frame=coord_frame)
    if is_j2000:
        result = result.transform_to(_fk5)
    return result

def get_albums():
    r = send_request("method_sync",{"method":"get_albums","params":""})
    return r.json()

def get_stack_setting():
    r = send_request("method_sync",{"method":"get_stack_setting"})
    return r.json()

def set_stack_setting(save_discrete_frame, save_discrete_ok_frame, light_duration_min=-1):
    r = send_request("method_sync",{"method":"set_stack_setting","params":{"save_discrete_frame":save_discrete_frame,"save_discrete_ok_frame":save_discrete_ok_frame,"light_duration_min":light_duration_min}})

def init(lat, lon, lang):
    send_request("method_sync",{"method":"pi_is_verified"})
    tz_name = tzlocal.get_localzone_name()
    tz = tzlocal.get_localzone()
    now = datetime.now(tz)
    date_json = {}
    date_json["year"] = now.year
    date_json["mon"] = now.month
    date_json["day"] = now.day
    date_json["hour"] = now.hour
    date_json["min"] = now.minute
    date_json["sec"] = now.second
    date_json["time_zone"] = tz_name
    date_data = {}
    date_data['method'] = 'pi_set_time'
    date_data['params'] = [date_json]
    send_request("method_sync",date_data)

    loc_data = {}
    loc_param = {}
    loc_param['lon'] = lon
    loc_param['lat'] = lat
    loc_param['force'] = False
    loc_data['method'] = 'set_user_location'
    loc_data['params'] = loc_param
    send_request("method_sync",loc_data)

    lang_data = {}
    lang_data['method'] = 'set_setting'
    lang_data['params'] = {'lang': lang}
    send_request("method_sync",lang_data)

def park():
    send_request("method_sync",{"method":"scope_park"})

def move_to_horizon():
    send_request("method_sync",{"method":"scope_move_to_horizon"})
