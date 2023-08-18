import requests
import json

configuration_data = {
    "long": 23.5,
    "lat": 10.0,
    "serre_surface1": 150.25
}

api_endpoint = 'http://127.0.0.1:5000/calculate_pv'

response = requests.post(api_endpoint, json=configuration_data)

if response.status_code == 200:

    result = response.json()
    print(result)
else:
    print("Error occurred. Status code:", response.status_code)
    print(response.text)
