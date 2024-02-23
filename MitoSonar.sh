#!/bin/bash

helpFunction()
{
   echo -e "***MitoSonar: Metabarcoding de Peixes***"
   echo -e ""
   echo -e "\tUso: $0 [options]"
   echo -e ""
   echo -e "FLAG\t\tOPT\t\tDESCRIPTION"
   echo -e "[Quick run]"
   echo -e "-d\t\tDefault\t\tRun with everything as default."
   echo -e ""
   echo -e "[Filtering parameters]"
   echo -e "-m\t\tmaxN\t\tDefault 0. After truncation, sequences with more than maxN Ns will be discarded. Note that dada currently does not allow Ns."
   echo -e "-q\t\ttruncQ\t\tDefault 2. Truncate reads at the first instance of a quality score less than or equal to truncQ."
   echo -e "-l\t\ttruncLen\tDefault 100. Truncate reads after truncLen bases. Reads shorter than this are discarded."
   echo -e "-t\t\ttrimLeft\tDefault 18. The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft."
   echo -e "-e\t\tmaxEE\t\tDefault 2. After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))"
   echo -e ""
   exit 1 # Exit script after printing help
}

maxN=0
truncQ=2
truncLen=100
trimLeft=18
maxEE=2

while getopts "hdm:ql:t:e:" opt; do
   case "$opt" in
      h ) helpFunction ;;
      d ) echo -e "[Running default settings]";;
      m ) maxN="${OPTARG:}" ;;
      q ) truncQ="${OPTARG:}" ;;
      l ) truncLen="${OPTARG:}" ;;
      t ) trimLeft="${OPTARG:}" ;;
      e ) maxEE="${OPTARG:}" ;;
      \? ) echo -e "Opção inválida: -$OPTARG" >&2; helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done



if [ $# -eq 0 ];
then
    helpFunction
fi

### MAIN SCRIPT

echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@                                                           @@@@@@@@@@"
echo "@@@@@@@@@                                                               @@@@@@@@"
echo "@@@@@@@@,                                  (###/                        @@@@@@@@"
echo "@@@@@@@@                                 #(     ###(#                   @@@@@@@@"
echo "@@@@@@@@                                   #####   .#((*                @@@@@@@@"
echo "@@@@@@@@                                        (((,  ///               @@@@@@@@"
echo "@@@@@@@@                          (%%%%@@    ##*  ///  ///              @@@@@@@@"
echo "@@@@@@@@                          %%%%%%%#@       .//  ,//              @@@@@@@@"
echo "@@@@@@@@                         #   .       . .@@@    .                @@@@@@@@"
echo "@@@@@@@@         %%%@&           *#*     #. #( ,%%%%&@@%                @@@@@@@@"
echo "@@@@@@@@          ##%%%%      (#     ## .(..(( ,%%%%%%%%@@              @@@@@@@@"
echo "@@@@@@@@           ,##(%%         ((.(( .(..(( .%#((   (#((#            @@@@@@@@"
echo "@@@@@@@@            .((## ./.  / .//.// ./..// .((((((((((/(            @@@@@@@@"
echo "@@@@@@@@          .(((//*     /(     //  /..//..(((////(/,              @@@@@@@@"
echo "@@@@@@@@         ((////          *#*    .*. // .//////                  @@@@@@@@"
echo "@@@@@@@@        ,,.                  //                                 @@@@@@@@"
echo "@@@@@@@@                           ///////.                             @@@@@@@@"
echo "@@@@@@@@                           .**.                                 @@@@@@@@"
echo "@@@@@@@@                                                                @@@@@@@@"
echo "@@@@@@@@      ( \/ )(  )(_  _)/  \   / ___) /  \ (  ( \ / _\ (  _ \     @@@@@@@@" 
echo "@@@@@@@@      / \/ \ )(   )( (  O )  \___ \(  O )/    //    \ )   /     @@@@@@@@"
echo "@@@@@@@@      \_)(_/(__) (__) \__/   (____/ \__/ \_)__)\_/\_/(__\_)     @@@@@@@@"
echo "@@@@@@@@,                                                               @@@@@@@@"
echo "@@@@@@@@@                                                               @@@@@@@@"
echo "@@@@@@@@@@@                                                           @@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"


## Chamada do script R:
Rscript MitoSonar_script.r "$maxN" "$truncQ" "$truncLen" "$trimLeft" "$maxEE"