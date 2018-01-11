import numpy as np

colloq_size=16
n_proc=4

for i in range(n_proc):
    #for colloq size
    start_index=i*(colloq_size/n_proc)
    if (i<(colloq_size % n_proc)):
        start_index=start_index+i
    else:
        start_index=start_index+(colloq_size % n_proc)
    end_index=start_index+(colloq_size/n_proc)-1
    if (i<(colloq_size % n_proc)):
        end_index=end_index+1
    if (i!=(n_proc-1)):
        end_index=end_index+1
    to_print="Colloq indices for processor "+str(i+1)+" is from "+str(start_index)+" to "+str(end_index)+". \n"
    print(to_print)

print ("\n")

for i in range(n_proc):
    #for stag size
    start_index=i*(colloq_size/n_proc)
    if (i<(colloq_size % n_proc)):
        start_index=start_index+i
    else:
        start_index=start_index+(colloq_size % n_proc)
    end_index=start_index+(colloq_size/n_proc)-1
    if (i<(colloq_size % n_proc)):
        end_index=end_index+1
    if (i!=(n_proc-1)):
        end_index=end_index+1
    if (i!=0):
        start_index=start_index-1
    end_index=end_index-1
    to_print="Staggered indices for processor "+str(i+1)+" is from "+str(start_index)+" to "+str(end_index)+". \n"
    print(to_print)
