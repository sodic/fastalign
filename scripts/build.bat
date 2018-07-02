@echo off
title builder

cmake -H. -Bbuild
cmake --build build