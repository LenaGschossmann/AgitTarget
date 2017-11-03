%EDFREAD Read Binary Eyelink Data Files
%
% TRIALS = EDFREAD(FILENAME)
%   Reads the EDF given by the string FILENAME.
%   Returns the Trials as entries in the TRIALS
%   Struct Array
%
% [TRIALS, INFO] = EDFREAD(FILENAME)
%   Additionally returns the Struct INFO with
%   Recording and Calibration Information
%
% 
% EDFREAD(FILENAME, FILTER*)
%   Every FILTER is a string to match for in the
%   message payloads in the EDF.
%   Every FILTER is added as a struct entry in
%   TRIALS.
%
% $Id: edfread.m 2 2007-06-17 11:09:53Z jsteger $
