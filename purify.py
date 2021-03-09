
# Extract code blocks from an .Rmd file

import sys, textwrap

print("# This file is generated from the corresponding .Rmd file")
print()
print()

in_code = False
in_challenge = False
counts = [ ]
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith("```"):
        print()
        in_code = not in_code
        assert in_code or line == "```", line
    elif in_code:
        print(line)
    elif line.startswith("#"):
        print("#" if in_challenge else "")

        in_challenge = "{.challenge}" in line

        if ".unlisted" in line: continue

        if in_challenge:
            line = line.replace("{.challenge}","").rstrip()
        #n = line.count("#")
        #line = line[n:].strip()
        #left = "#"*n+" "
        #bracket = "----" if n > 1 else "===="
        #print(left+"_"*(n-1+len(line)+len(bracket)*2+2))
        #print(left+bracket+">"*(n-1)+" "+line+" "+bracket)

        n = line.count("#")
        while n < len(counts): del counts[-1]
        while n > len(counts): counts.append(0)
        counts[-1] += 1

        banner = "# "+".".join(str(i) for i in counts) + line[n:] + " ----"
        if n == 1:
            print()
            print()
            print("#"+"/"*(len(banner)-1))
        print(banner)

    elif in_challenge:
        for line2 in textwrap.wrap(line) or [""]:
            print("# " + line2)
