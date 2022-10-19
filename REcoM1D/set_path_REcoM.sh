#!/bin/bash
cd .. 
dir=$(pwd)
export RECOM_MAIN_DIR=$dir
export RECOM_NAMELIST_PATH=${RECOM_MAIN_DIR}/config/
export RECOM_GRID_PATH=${RECOM_MAIN_DIR}/grid/
export RECOM_DATA_PATH=${RECOM_MAIN_DIR}/data/
export RECOM_FORCING_PATH=${RECOM_MAIN_DIR}/forcing/
export RECOM_RESULT_PATH=${RECOM_MAIN_DIR}/results/
cd bin

