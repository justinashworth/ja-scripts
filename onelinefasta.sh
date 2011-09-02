#!/bin/bash
cat $1 | sed -e 's/>\(.*\)/>\1#/g' | tr -d '\n' | tr '>' '\n>' | tr '#' '\n'
