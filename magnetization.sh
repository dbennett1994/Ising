#!/bin/bash

echo "Please enter values of J in range (0,5] separated by a space."

read -a values
echo "Thank you. Please wait for the graph of magnetization vs temperature to be produced."

./magnetization.py echo ${values[@]}

