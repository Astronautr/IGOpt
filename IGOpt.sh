#!/bin/bash

javac -classpath ./lib/p1788-rev444.jar -sourcepath ./src -d ./out ./src/Main.java
java -classpath ./lib/p1788-rev444.jar:./out Main $@
