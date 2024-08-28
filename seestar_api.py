import requests
import json
import time
from astropy.coordinates import SkyCoord

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

def move_ra():
    total_move = 90
    subtotal_move = 0
    while subtotal_move < total_move:
        cur_move = min(45, total_move - subtotal_move)
        move_time = cur_move * 20.0 / 90.0
        send_request("method_sync", {"method":"scope_speed_move","params":{"speed":1000,"angle":0,"dur_sec":round(move_time)}})
        subtotal_move += cur_move
        time.sleep(move_time + 1)

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
    # don't know why there is a \ in the path
    send_request("method_sync",{"method":"delete_image", "params":{"name":[f"MyWorks\\/{path}"]}})

def goto(target_name, ra, dec, is_j2000=True):
    # ra and dec are in hms and dms respectively
    parameters = {"target_name":target_name, "ra":ra, "dec":dec, "is_j2000":is_j2000}
    send_request("goto_target",parameters)

def goto_without_platesolve(target_name, ra, dec, is_j2000=True):
    # ra and dec are in decimal hour and degree respectively
    set_target_name(target_name)
    target_coords = SkyCoord(ra=ra, dec=dec)
    parameters = {"method":"scope_goto","params":[target_coords.ra.hour,target_coords.dec.degree]}
    send_request("method_sync",parameters)

def wait_for_goto():
    print("Waiting for goto to complete")
    while is_goto():
        time.sleep(2)

def wait_for_goto_ra_dec(target_ra, target_dec):
    print("Waiting for goto ra and dec")
    coords = get_equ_coords()
    while (coords["ra"] - target_ra) ** 2 + (coords["dec"] - target_dec) ** 2 > 0.01:
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
