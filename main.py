from alignment import *

A = input("Enter first string: ")
B = input("Enter second string: ")

ali = Alignment(-5)
ali.local_alignment(A, B)
