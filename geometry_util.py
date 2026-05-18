import json
import re


geo_to_IC0 = 1948.07 #distance from surface of IC geometry to IC z=0

device_name_map = {"mdom":"mDOM","degg":"DEgg","pdom":"pDOM","lom":"LOM"}

def get_device_type(icm_id,geometry_list):
    device_type = None
    for geometry_file in geometry_list:
        with open(geometry_file, 'r') as f:
            geometry_data = json.load(f)
            for device in geometry_data[0]['devices']:
                if device['icm_id'] == icm_id:
                    device_type = device_name_map[device['type'].lower()]
    return device_type


def get_string_port_from_icm_id(icm_id,string_map_list):
    string = None
    port = None
    for string_map in string_map_list:
        with open(string_map, 'r') as f:
            data = json.load(f)
            for device in data:
                # print(device)
                if device["icm_id"] == icm_id:
                    hostname = device["hostname"]
                    port = device["port"]
                    string = re.search(r'string_(\d+)_', string_map).group(1)
    return string, port

def get_depth_from_icm_id(icm_id,geometry_list):
    depth = None
    for geometry_file in geometry_list:
        with open(geometry_file, 'r') as f:
            geometry_data = json.load(f)
            for device in geometry_data[0]['devices']:
                if device['icm_id'] == icm_id:
                    depth = geo_to_IC0 - device['z']
    return depth

def get_string_position_from_icm_id(icm_id,geometry_list):
    position = None
    for geometry_file in geometry_list:
        with open(geometry_file, 'r') as f:
            geometry_data = json.load(f)
            for device in geometry_data[0]['devices']:
                if device['icm_id'] == icm_id:
                    position = device['position']
                    upgrade_string = geometry_data[0]['id']
    return position,upgrade_string


def prod_id_to_icm_id(prod_id,geometry_files) -> str:
    '''
    production id
    '''
    icm_id_list = []
    for ifile in geometry_files:
        with open(ifile, 'r') as f:
            data = json.load(f)
        for idev in data[0]["devices"]:
            if idev["production_id"] == prod_id:
                icm_id_list.append(idev["icm_id"])
    if len(icm_id_list)>1:
        print(f"more than one icm id found for prod id {prod_id} {icm_id_list}")
    return icm_id_list[0]


def get_icm_id_from_string_port(string,port,string_map_list):
    string_map_file = [istring_map for istring_map in string_map_list if f"string_{string}_" in istring_map][0]
    icm_id = None
    with open(string_map_file, 'r') as f:
        data = json.load(f)
        for device in data:
            if device["hostname"] == f"fieldhub{string}" and device["port"] == str(port):
                icm_id = device["devices"][0]["icm_id"]
    return icm_id

def get_string_device_from_string_port(string,port,string_map_list,geometry_list):
    return get_string_position_from_icm_id(get_icm_id_from_string_port(string, port, string_map_list), geometry_list)