import requests
import json
import time

def send_request(action, parameters):
    ip_address = "localhost"
    port = 5555
    device_num = 1
    base_url = f"http://{ip_address}:{port}/api/v1/telescope/{device_num}/action"

    headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"}
    body = {"Action":action,"Parameters":json.dumps(parameters),"ClientID":1,"ClientTransactionID":999}
    r = requests.put(base_url,body,headers=headers)
    print(body)
    print(r.status_code)
    print(r.text)

def move_ra():
        total_move = 90
        subtotal_move = 0
        #        self.sync_target([0.0, -89.0])
        while subtotal_move < total_move:
            cur_move = min(45, total_move - subtotal_move)  # 45 = 90/20/10
            move_time = cur_move * 20.0 / 90.0
            send_request("method_sync", {"method":"scope_speed_move","params":{"speed":1000,"angle":180,"dur_sec":round(move_time)}})
            #self.move_scope(0, 1000, round(move_time))
            subtotal_move += cur_move
            time.sleep(move_time + 1)
            #       self.sync_target([0.0, 80.0])

move_ra()