import matplotlib.pyplot as plt
files = ['4', '16', '64', '128']

with open('1', 'r') as f:
	serial_time = float(f.read().split("time: ")[-1].replace('\t', '')) 
points = [(1, 1.0)]
for name in files:
	total_time = 0
	count = 0
	with open(name, 'r') as f:
		for line in f:
			if 'time: ' in line:
				total_time += float(line.split('time: ')[-1].replace('\t', ''))
				count += 1
	points.append((int(name), serial_time / (total_time / count)))

print(points)

