#!/bin/bash
for line in `cat $1`;do
  ssh $line "echo $line; ps ax |pgrep java"
done

