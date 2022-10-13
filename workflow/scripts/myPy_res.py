#!/usr/bin/env python
# -*- coding: utf-8
 

## set time resource of a rule - VSC5
def time_MEGAHIT(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 4:
        time = ("00:%02d:00" % (59))
    elif in_size < 20:
        time = ("%02d:00:00" % (12))
    elif in_size < 45:
        time = ("%02d:00:00" % (24))
    else:
        time = ("%02d:00:00" % (32))

    return time

def qos_MEGAHIT(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 45:
        out_qos = "p71863_2048"
    else:
        out_qos = "p71863_0512"

    return out_qos

def part_MEGAHIT(wildcards, input):    
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 45:
        out_part = "zen3_2048"
    else:
        out_part = "zen3_0512"
    
    return out_part


# def node_HIPMER(wildcards, input):
#     in_size = input.size_mb/1024
#     in_size = max(1, int(in_size))

#     if in_size < 15:
#         node = 1
#     elif in_size < 40:
#         node = 2
#     elif in_size < 120:
#         node = 1
#     else:
#         node = 2

#     return node

def node_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 20:
        node = 1
    elif in_size < 40:
        node = 2
    elif in_size < 100:
        node = 3
    elif in_size < 140:
        node = 4
    elif in_size < 170:
        node = 5
    else:
        node = 6

    return node

def time_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 15:
        time = ("00:%02d:00" % (59))
    elif in_size < 40:
        time = ("00:%02d:00" % (59))
    elif in_size < 120:
        time = ("%02d:00:00" % (4))
    else:
        time = ("%02d:00:00" % (12))

    return time

def qos_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size < 15:
        qos = "p71863_0512"
    elif in_size < 40:
        qos = "p71863_0512"
    elif in_size < 120:
        qos = "p71863_2048"
    else:
        qos = "p71863_2048"

    return qos

def part_HIPMER(wildcards, input):
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))
    
    if in_size < 15:
        part = "zen3_0512"
    elif in_size < 40:
        part = "zen3_0512"
    elif in_size < 120:
        part = "zen3_2048"
    else:
        part = "zen3_2048"

    return part
## ---










### deprecated:


## set time resource of a rule - VSC4
def get_time_spades(wildcards, input):
    import math

    time_coef = 7
    run_time = input.size_mb/1024 * time_coef
    hours = max(2, math.ceil(run_time)) # to ensure min time of 2h
    if hours > 72:
        return -1 # kill job if it takes more than 3d
    time = ("%02d:00:00" % (hours))
    return time

def get_time_megahit(wildcards, input):
    import math
    
    time_coef = 0.5
    run_time = input.size_mb/1024 * time_coef
    hours = max(1, math.ceil(run_time)+1) # to ensure min time of 1h
    time = ("%02d:00:00" % (hours))
    return time

def get_qos_vsc4(wildcards, input): # VSC4
    in_size = input.size_mb/1024
    in_size = max(1, int(in_size))

    if in_size > 30:
        out_qos = "mem_0384"
    else:
        out_qos = "mem_0096"

    return out_qos
## ---