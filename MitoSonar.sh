#!/bin/bash

helpFunction()
{
   echo -e "***MitoSonar: Multipurpose Metabarcoding Tool***"
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
echo "@@@@@@@@       _  _  __  ____  __     ____   __   __ _   __   ____      @@@@@@@@"
echo "@@@@@@@@      ( \/ )(  )(_  _)/  \   / ___) /  \ (  ( \ / _\ (  _ \     @@@@@@@@" 
echo "@@@@@@@@      / \/ \ )(   )( (  O )  \___ \(  O )/    //    \ )   /     @@@@@@@@"
echo "@@@@@@@@      \_)(_/(__) (__) \__/   (____/ \__/ \_)__)\_/\_/(__\_)     @@@@@@@@"
echo "@@@@@@@@,                                                               @@@@@@@@"
echo "@@@@@@@@@                                                               @@@@@@@@"
echo "@@@@@@@@@@@                                                           @@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"


## Interface de seleção de banco de dados
while true; do
  echo "Please, select the kind of tax-search you wish to perform:"
  echo "   1) Fish mitogenomes (12s-like)"
  echo "   2) Bacteria ribossomal RNA (16s-like)"
  echo "   3) Fungi ribossomal RNA (ITS-like)"
  read -p "Enter your choice: " choice

  case $choice in
    1)
      echo "You selected Fish mitogenomes (12s-like)"
      database='MiFish'
      break
      ;;
    2)
      echo "You selected Bacteria ribossomal RNA (16s-like)"
      database='Silva'
      break
      ;;
    3)
      echo "You selected Fungi ribossomal RNA (ITS-like)"
      database='UNITE'
      break
      ;;
    *)
      printf "\nWARNING: Invalid selection. Please enter the number of your choice.\n\n"
      ;;
  esac
done

## Definindo se amostra é de fita simples ou dupla
while true; do
  echo "Please, inform if your input data is single or pair-ended"
  echo "   1) Single strand"
  echo "   2) Pair-ended"
  read -p "Enter your choice: " choice

  case $choice in
    1)
      echo "Initiating tool for simple input..."
      inputtype='single'
      break
      ;;
    2)
      echo "Initiating tool for pair-ended input..."
      inputtype='pair'
      break
      ;;
    *)
      printf "\nWARNING: Invalid selection. Please enter the number of your choice.\n\n"
      ;;
  esac
done


## Chamada do script R:
Rscript MitoSonar_script.r "$maxN" "$truncQ" "$truncLen" "$trimLeft" "$maxEE" "$database" "$inputtype"