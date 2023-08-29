import requests
import json

configuration_data = {
    "long": -1.0285,
    "lat": 43.5161,
}

api_endpoint = 'http://api.enervivo.fr/calculate_pv'

response = requests.post(api_endpoint, json=configuration_data)

if response.status_code == 200:

    result = response.json()
    print(result)
else:
    print("Error occurred. Status code:", response.status_code)
    print(response.text)
