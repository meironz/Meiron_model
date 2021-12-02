with open('environment.yml','r') as f:
    lines = f.readlines()

newlines = []    
for line in lines:
    newlines.append(line.split('=')[0]+'\n')
    
with open('clean_environment.yml','w') as f:
    f.writelines(newlines)


