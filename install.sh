#!/usr/bin/env bash

virtualenv -p python2.7 --no-site-packages --always-copy env
source env/bin/activate
pip install --no-cache-dir --ignore-installed --force-reinstall  .


