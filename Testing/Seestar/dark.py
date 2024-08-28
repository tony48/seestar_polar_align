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

r = send_request("method_sync",{"method":"set_setting","params":{"frame_calib":True}})
r = send_request("method_sync",{"method":"set_setting","params":{"auto_3ppa_calib":True}})