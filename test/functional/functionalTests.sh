#!/bin/bash

host=`hostname -s`

if [ "$host" == "lam" ]; then
    ./linematching2.py
    exit 0
fi

./linematching2.py > /home/jenkins/jenkins/amazed/cpf-redshift/label/linux-centos7/test/report/functional.xml
