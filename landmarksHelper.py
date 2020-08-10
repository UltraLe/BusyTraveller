def readCSV(filename):
	fd = open(filename, "r")

	lines = fd.readlines()

	names = lines[0].strip().split(", ")

	print(names)
	final_list = []
	for i in range(1, len(lines)):
		values = lines[i].strip().split(", ")
		if len(values) != len(names):
			values = lines[i].strip().split(",")
		dict_info = {}	
		for j in range(len(names)):
			dict_info[names[j]] = values[j].strip()
		final_list.append(dict_info)

	return final_list


def recover_landmarks():

	#@ return: list of monuments names

	res = []

	fd = open("landmarks.csv", "r")

	lines = fd.readlines()
	lines = lines[1:]
	for line in lines:
		splitted = line.strip().split(", ")
		res.append(splitted[0]) 
	fd.close()

	return res[:50]

def recover_distances():

	#@ return: list of dict with key "Start" "End" "Seconds" and "Transport"

	res = []

	fd = open("distance.csv", "r")

	lines = fd.readlines()

	for i in range(1, len(lines)):
	    splitted = lines[i].strip().split(", ")
	    res.append({"Start": splitted[0], "End": splitted[1], "Seconds": splitted[2], "Transport": splitted[3]})
	fd.close()

	return res
