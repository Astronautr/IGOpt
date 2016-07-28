#!/bin/bash

# if [[ $# -gt 0 ]]; then
	
# 	files=$(find . ~ -maxdepth 1 -name $1)
# 	files=${files//"./"}
# 	if [[ "$files" == "$1" ]]; then
# 		echo "input file found"
# 	else
# 		echo "error: input file not found"
# 	fi

# else
# 	echo "Conditions will be grab from input"
# fi

# if [[ "$str1" == "$str2" ]]; then
# 	echo "equals"
# fi
javac -classpath ./lib/p1788-rev444.jar -sourcepath ./src -d ./out ./src/Main.java
java -classpath ./lib/p1788-rev444.jar:./out Main $@
