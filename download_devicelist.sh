

export BASE_URL=https://verical.icecube.wisc.edu

VERICAL_USER="icecube"
VERICAL_PASS="skua"

# curl -s -u $VERICAL_USER:$VERICAL_PASS -o monidaq_fieldhub_fieldhub87.json "$BASE_URL/monidaq/fieldhub/fieldhub87/raw?download=1"
# curl -s -u $VERICAL_USER:$VERICAL_PASS "$BASE_URL/monidaq/device-list/raw" | jq '.[].devices[].icm_id'
curl -s -u $VERICAL_USER:$VERICAL_PASS -o monidaq_device_list.json "$BASE_URL/monidaq/device-list/raw?download=1"