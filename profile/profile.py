from math import log

msa = {}
length = 0

n = int(input())
for i in range(n):
    sequence = input()
    chars = list(sequence)
    length = len(chars)
    for j, ch in enumerate(chars):
        if msa.get(ch) is None:
            msa[ch] = [j]
        else:
            positions = msa.get(ch)
            positions.append(j)

total = 2 * len(msa) * length
probabilities = {}
for key in msa:
    counts = [1] * length
    sum = length + len(msa.get(key))
    for position in msa.get(key):
        counts[position] += 1
    quotients = [count / sum for count in counts]
    logs = [log(y, 10) for y in quotients]
    probabilities[key] = logs

req = input()
chars = list(req)
max_score = -10000
max_string = []
for i in range(len(chars)-length):
    sub = chars[i:i+length]
    score = 0
    for j, ch in enumerate(sub):
        score += probabilities.get(ch)[j]
    if score > max_score:
        max_score = score
        max_string = sub

ans = ""
print(ans.join(max_string))
