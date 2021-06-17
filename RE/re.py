expression = input().split("-")
n = int(input())
for i in range(n):
    seq = input()
    counter = 0
    error = 0
    match_string = ""
    for p in expression:
        repeat = 1
        if p[len(p) - 1] == ")":
            repeat = int(p[len(p)-2:len(p)-1])
            p = p[:len(p)-2]
        for _ in range(repeat):
            # a letter
            if len(p) == 1:
                if seq[counter] == p:
                    match_string += "1"
                else:
                    error += 1
                    match_string += "0"
                counter += 1
            elif p[0] == "[":
                letters = p[1:len(p)-1]
                match = False
                for l in letters:
                    if l == seq[counter]:
                        match_string += "1"
                        match = True
                        break
                counter += 1
                if not match:
                    match_string += "0"
                    error += 1
            else:
                letters = p[1:len(p) - 1]
                match = False
                for l in letters:
                    if l == seq[counter]:
                        match_string += "0"
                        error += 1
                        match = True
                        break
                counter += 1
                if not match:
                    match_string += "1"
    if error > 2:
        print("No ", match_string)
    else:
        print("Yes ", match_string)
