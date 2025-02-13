#!/bin/bash

ml idr/2.0.3

idr --samples CQTL_AM7754_CTL_Ankle_1_peaks.bed CQTL_AM7754_CTL_Ankle_replicate_peaks.bed \
	--input-file-type bed \
	--rank 
