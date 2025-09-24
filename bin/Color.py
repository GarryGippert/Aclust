''' package Color.py written by Garry P. Gippert '''

import re, colorsys, struct, binascii

morecolors = [
        'F2460C',
        '6EAE83',
        'BC9EF9',
        '3D3D8E',
        '653D8E',
        '8E8E3D',
        '0CF2B8',
        'F9BC9E',
        'F99EDB',
        '7F0CF2',
        '0CF246',
        'B80CF2',
        '46F20C',
        'F20CB8',
        'F20C0C',
        '9EF9F9',
        'F2F20C',
        '0C46F2',
        '9EBCF9',
        '668E3D',
        '8E3D66',
        '9EF9DB',
        'B8F20C',
        '3D8E8E',
        'BCF99E',
        '9E9EF9',
        'F2B80C',
        'F99EBC',
        'F20C7F',
        '8E663D',
        'DBF99E',
        '460CF2',
        'F9F99E',
        '9EF99E',
        '8E3D8E',
        '3D8E66',
        'F99E9E',
        '0CF27F',
        'F99EF9',
        'F20C46',
        '9EDBF9',
        '0C7FF2',
        '0CB8F2',
        'F9DB9E',
        '0C0CF2',
        'F20CF2',
        '0CF2F2',
        '7FF20C',
        '3D8E3D',
        '0CF20C',
        '3D658E',
        '8E3D3D',
        'DB9EF9',
        'F27F0C', ]

def _hex2rgb(hex):
	hex = hex.upper()
	hex = re.sub('^#', '', hex)
	rgb = struct.unpack('3B', binascii.unhexlify(hex))
	return rgb

def _rgb2hex(rgb):
	hex = "%02X%02X%02X" % ( rgb[0], rgb[1], rgb[2] )
	return hex

def _hexbg2fg(hex):
	# return black 000000 or white ffffff based on best contrast to background color
	rgb = _hex2rgb(hex)
	y = 0.2126 * rgb[0] / 256 + 0.7152 * rgb[1] / 256 + 0.0722 * rgb[2] / 256;
	return "000000" if y > 0.5  else "FFFFFF"

def palecolor(hex, f=0.5):
	''' add fraction f of lightness to a hex color code and return '''
	if abs(f)>1.0:
		raise Exception("Color:palecolor", hex, 'fraction abs(', f, ') > 1.0')
	# convert hex to RGB
	R, G, B = _hex2rgb(hex)
	# convert RGB to HLS
	H, L, S = colorsys.rgb_to_hls(R/255.0, G/255.0, B/255.0)
	# increase or decrease lightness by f * (1-L)
	L = max(L + f * (1.0 - L), 0.0)
	# convert new HLS to RGB
	R, G, B = colorsys.hls_to_rgb(H, L, S)
	rgb = [int(R*255.0), int(G*255.0), int(B*255.0)]
	# convert new RGB to hex
	return _rgb2hex(rgb)

def zeropad(str, pad=10):
	''' general zero padding of strings containing substrings integers
		for example: GH15a_17 -> GH0000000015a_0000000017 '''
	w = re.split('([0-9]+)', str)
	fmt = f"%0{pad}d"
	for i in range(0, len(w)):
		try:
			w[i] =  fmt % int(w[i])
		except:
			pass
	return ''.join(w)
