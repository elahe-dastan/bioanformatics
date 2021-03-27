from alignment import *

A = input("Enter first string: ")
B = input("Enter second string: ")

ali = Alignment(-5, A, B)
ali.local_alignment()
