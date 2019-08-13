#!/bin/bash

declare -a cols  ## array holding original columns from original data file
declare -a csel  ## array holding columns to select (from file 2)
declare -a cpos  ## array holding array indexes of matching columns

cols=( $(head -n 1 "$1") )  ## fill cols from 1st line of data file
csel=( $(< "$2") )          ## read select columns from file 2

## fill column position array
for ((i = 0; i < ${#csel[@]}; i++)); do
    for ((j = 0; j < ${#cols[@]}; j++)); do
        [ "${csel[i]}" = "${cols[j]}" ] && cpos+=( $j )
    done
done

printf " " 
for ((i = 0; i < ${#csel[@]}; i++)); do   ## output header row
    printf "    %s" "${csel[i]}"
done

printf "\n"     ## output newline
unset cols      ## unset cols to reuse in reading lines below

while read -r line; do        ## read each data line in data file 
    cols=( $line )            ## separate into cols array
    printf "%s" "${cols[0]}"  ## output row label
    for ((j = 0; j < ${#cpos[@]}; j++)); do
        [ "$j" -eq "0" ] && { ## handle format for first column
            printf "%5s" "${cols[$((${cpos[j]}))]}"
            continue
        }                     ## output remaining columns
        printf "%13s" "${cols[$((${cpos[j]}))]}"
    done
    printf "\n"
done < <( tail -n+2 "$1" )

