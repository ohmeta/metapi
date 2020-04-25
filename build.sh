#!/bin/sh
rm -rf build
rm -rf metapi.egg-info
rm -rf dist
python3 setup.py sdist bdist_wheel
twine upload dist/*
