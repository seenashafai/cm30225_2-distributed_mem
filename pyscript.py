import subprocess

n = [10, 20]

arr = []
fout = open("./Output.txt", "w")


for r in range(0, len(n)):


    for i in range(0, n[r]*n[r]):
      if (i < n[r] or ((i) % n[r] == 0) or ((i+1) % n[r] == 0) or (i>= (n[r]*n[r])-n[r])):
          arr.insert(i, 1)
      else:
          arr.insert(i, 0) #Otherwise assign 0


    with open('textfile.txt', 'w') as f:
        for item in arr:
            f.write("%s\n" % item)

    subprocess.Popen(['mpirun', 'dist_mem', str(n[r]), '0.01'], stdout=fout)
    arr = []


fout.close()
