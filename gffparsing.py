import sys,os
from BCBio import GFF
from Bio import SeqIO

inputtype=sys.argv[1]
#outpath=sys.argv[2]
#out_handle = open(outpath, "w")

if inputtype=='file':
    in_file=sys.argv[2]
    #write gff
    in_handle = open(in_file)
    GFF.write(SeqIO.parse(in_handle, "genbank"), sys.stdout)
    in_handle.close()
    #out_handle.close()
else:
    GFF.write(SeqIO.parse(sys.stdin, "genbank"), sys.stdout)
    #out_handle.close()




#OLD CODE

#genomename=sys.argv[3]

##check file has features
#with open(in_file) as f:
#    for indx, seq_record in enumerate(SeqIO.parse(f, "genbank")):
#        if seq_record.features:
#            pass
#        else:
#            print('Error: genbank file has no feature annotation - are you using the Genbank full format?')
#            sys.exit()


# #check path is valid file #do this check in driver script
# if (os.path.isfile(in_file)):
#     pass
# else:
#     print('Error: %s is not a file'%in_file)
#     sys.exit()
