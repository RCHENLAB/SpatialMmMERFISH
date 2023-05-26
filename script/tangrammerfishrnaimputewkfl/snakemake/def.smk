import sys

## Type convertion from string to others
def opvalue(op, value):
	if op in ['=']:
		return value
	elif op in ['s', 'i', 'f', 'b']:
		if op=='s':
			return value
		elif op=='i':
			return int(value)
		elif op=='f':
			return float(value)
		elif op=='b':
			return bool(value)
	elif op in ['S', 'I', 'F', 'B']:
		if op=='S':
			return value.split('::')
		elif op=='I':
			return [int(v) for v in value.split('::')]
		elif op=='F':
			return [float(v) for v in value.split('::')]
		elif op=='B':
			return [bool(v) for v in value.split('::')]
	else:
		print(f"Error: {op} is not recognized.")
		sys.exit(-1)
