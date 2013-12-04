oneletter.key = list(
	ALA = 'A',
	CYS = 'C',
	ASP = 'D',
	GLU = 'E',
	PHE = 'F',
	GLY = 'G',
	HIS = 'H',
	ILE = 'I',
	LYS = 'K',
	LEU = 'L',
	MET = 'M',
	ASN = 'N',
	PRO = 'P',
	GLN = 'Q',
	ARG = 'R',
	SER = 'S',
	THR = 'T',
	VAL = 'V',
	TRP = 'W',
	TYR = 'Y',

	A   = 'A',
	C   = 'C',
	D   = 'D',
	E   = 'E',
	F   = 'F',
	G   = 'G',
	H   = 'H',
	I   = 'I',
	K   = 'K',
	L   = 'L',
	M   = 'M',
	N   = 'N',
	P   = 'P',
	Q   = 'Q',
	R   = 'R',
	S   = 'S',
	T   = 'T',
	V   = 'V',
	W   = 'W',
	Y   = 'Y'
)

threeletter.key = list(
	A   = 'ALA',
	C   = 'CYS',
	D   = 'ASP',
	E   = 'GLU',
	F   = 'PHE',
	G   = 'GLY',
	H   = 'HIS',
	I   = 'ILE',
	K   = 'LYS',
	L   = 'LEU',
	M   = 'MET',
	N   = 'ASN',
	P   = 'PRO',
	Q   = 'GLN',
	R   = 'ARG',
	S   = 'SER',
	T   = 'THR',
	V   = 'VAL',
	W   = 'TRP',
	Y   = 'TYR',

	ALA = 'ALA',
	CYS = 'CYS',
	ASP = 'ASP',
	GLU = 'GLU',
	PHE = 'PHE',
	GLY = 'GLY',
	HIS = 'HIS',
	ILE = 'ILE',
	LYS = 'LYS',
	LEU = 'LEU',
	MET = 'MET',
	ASN = 'ASN',
	PRO = 'PRO',
	GLN = 'GLN',
	ARG = 'ARG',
	SER = 'SER',
	THR = 'THR',
	VAL = 'VAL',
	TRP = 'TRP',
	TYR = 'TYR'
)

oneletter =
	function(x)
{
	return( as.character( sapply(x, function(y) { oneletter.key[y] } ) ) )
}

threeletter =
	function(x)
{
	return( as.character( sapply(x, function(y) { threeletter.key[y] } ) ) )
}
