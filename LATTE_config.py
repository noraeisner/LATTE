
import os
import csv

# check what the current path is - when the program is first downloaded the path is set to 'no/path/set/yet' and the user is automatically prompted to change'no/path/set/yet'
with open("_config.txt", 'r') as f:
    path = str(f.readlines()[-1])

print (path)
new_path = True

def yes_or_no():

    '''
    Yes/No command line option to verify that the user really wants to change the output/input path
    '''
    print ('\n \n WARNING: if you have already downloded the input files (with --new-data) then these will remain in the location set by your previous path, so you will have to redowload the data (not recommended) or move the data to the new location set by this path. \n \n ')

    reply = str(input('Are you sure that you want to change the path?' + '(yes/no): '))

    if (reply == 'y') or (reply == 'yes') or (reply == 'yep') or (reply == 'yeah'):
        return True

    else: # if anything else is entered assume that this is a 'no' and continue with the old path
        return False
        


if path == 'no/path/set/yet':
    indir = input("\n \n No output path has been set yet. \n \n Please enter a path to save the files (e.g. ./LATTE_output or /Users/yourname/Desktop/LATTE_output) : " )

    # SAVE the new output path
    with open("_config.txt",'w') as f:
        f.write(str(indir))


    print("\n New path: " + indir)



# if the user chooses to redefine the path

elif new_path == True: 

    reply = yes_or_no()

    if reply == True:
        indir = input("\n \n Please enter a path to save the files (e.g. ./LATTE_output or /Users/yourname/Desktop/LATTE_output) : " )

        # SAVE the new output path
        with open("_config.txt",'w') as f:
            f.write(str(indir))    
        
        print("\n New path: " + indir)

    else:
        

        print ("LATTE will continue to run with the old path: {}".format(path))
        indir = path
else:
    indir = path











