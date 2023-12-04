#!/bin/sh

help="
  Usage: $(basename "$0") --i <invcf> --g <ingff> --f <infasta> --o <outpath>
  
  where:
     -i | --invcf   Specify the input VCF file
     -g | --ingff   Specify the input GFF file
     -f | --infasta Specify the input FASTA file
     -o | --outpath Specify the output path
     -h | --help    Show this help message
"
terminal=$(stty -g)

SHORT=i:g:f:o:h
LONG=invcf:,ingff:,infasta:,outpath:,help
OPTS=$(getopt -o $SHORT --long $LONG --name "$0" -- "$@")
if [ $? -ne 0 ]; then echo "$help" >&2; exit 1; fi
eval set -- "$OPTS"

while true; do
 case "$1" in
   -h | --help )    echo "$help"; exit 0  ;;
   -i | --invcf )   invcf="$2"  ; shift 2 ;;
   -g | --ingff )   ingff="$2"  ; shift 2 ;;
   -f | --infasta ) infasta="$2"; shift 2 ;;
   -o | --outpath ) outpath="$2"; shift 2 ;;
   -- )             shift; break ;;
 esac
done

for file in "$invcf" "$ingff" "$infasta" "$outpath"; do
  if [ -z "$file" ]; then
     echo "Error: There are some missing or invalid arguments." >&2
     echo "$help" >&2; exit 1
  fi
done

for file in "$invcf" "$ingff" "$infasta"; do
  if [ ! -f "$file" ]; then
     echo "Error: The file $(realpath "$file") does not exist." >&2; exit 1
  fi
done

mkdir -p "$outpath"
if [ ! -d "$outpath" ]; then
  echo "Error: The directory $(realpath "$outpath") does not exist." >&2; exit 1
fi

rmkdir() {
  if [ -d "$1" ] && [ ! -z "$(ls "$1")" ]; then
    echo "Warning: The directory $(realpath "$1") already exists and will be overwritten."
    echo -n "Do you want to proceed? (y/n)"
    while true; do
      stty raw -echo
      yn=$(dd bs=1 count=1 2>/dev/null)
      stty "$terminal"
      case "$yn" in 
        y|Y ) echo; break    ;;
        n|N ) echo; return 1 ;;
      esac
    done
    rm -rf "$1"/*
  fi
  mkdir -p "$1"
}

rmkdirs() {
  if [ -d "$1" ] && [ ! -z "$(ls "$1")" ] || [ -d "$2" ] && [ ! -z "$(ls "$2")" ]; then
    echo "Warning: The directories $(realpath "$1") and $(realpath "$2") already exist and will be overwritten."
    echo -n "Do you want to proceed? (y/n)"
    while true; do
      stty raw -echo
      yn=$(dd bs=1 count=1 2>/dev/null)
      stty "$terminal"
      case "$yn" in 
        y|Y ) echo; break    ;;
        n|N ) echo; return 1 ;;
      esac
    done
    rm -rf "$1"/* "$2"/*
  fi
  mkdir -p "$1" "$2"
}

if rmkdir "$outpath/VCF_Preprocessing"; then
  if [ $(awk -F'\t' 'NR==1{print NF}' "$invcf") -lt 8 ]; then
    awk 'BEGIN{ FS = OFS = "\t" } { print $0, "." }' "$invcf" > "$outpath/VCF_Preprocessing/$(basename "$invcf")"
  else
    cp "$invcf" "$outpath/VCF_Preprocessing/$(basename "$invcf")"
  fi
  python3 "$(dirname "$0")/navip.py" --mode pre --invcf "$outpath/VCF_Preprocessing/$(basename "$invcf")" --outpath "$outpath/VCF_Preprocessing/"
  rm "$outpath/VCF_Preprocessing/$(basename "$invcf")"
fi

if [ -f "$outpath/VCF_Preprocessing/first.vcf" ] && [ -f "$outpath/VCF_Preprocessing/second.vcf" ]; then
  # test if the preprocessing step actually produced two different VCF files for further processing
  if cmp -s "$outpath/VCF_Preprocessing/first.vcf" "$outpath/VCF_Preprocessing/second.vcf"; then
    
    if rmkdir "$outpath/NAVIP_Main_Output"; then
      rm -rf "$outpath/NAVIP_First_Output" "$outpath/NAVIP_Second_Output"
      python3 "$(dirname "$0")/navip.py" --mode main --invcf "$outpath/VCF_Preprocessing/first.vcf" --ingff "$ingff" --infasta "$infasta" --outpath "$outpath/NAVIP_Main_Output/"
    fi
    if [ -f "$outpath/NAVIP_Main_Output/All_VCF.vcf" ]; then
      if rmkdir "$outpath/SFA_Output"; then
        rm -rf "$outpath/SFA_First_Output" "$outpath/SFA_Second_Output"
        python3 "$(dirname "$0")/navip.py" --mode sfa --innavipvcf "$outpath/NAVIP_Main_Output/All_VCF.vcf" --outpath "$outpath/SFA_Output/"
      fi
    fi
    
  else # if the two VCF files are different, run the NAVIP main program for both files separately
    
    if rmkdirs "$outpath/NAVIP_First_Output" "$outpath/NAVIP_Second_Output"; then
      rm -rf "$outpath/NAVIP_Main_Output"
      python3 "$(dirname "$0")/navip.py" --mode main --invcf "$outpath/VCF_Preprocessing/first.vcf" --ingff "$ingff" --infasta "$infasta" --outpath "$outpath/NAVIP_First_Output/"
      python3 "$(dirname "$0")/navip.py" --mode main --invcf "$outpath/VCF_Preprocessing/second.vcf" --ingff "$ingff" --infasta "$infasta" --outpath "$outpath/NAVIP_Second_Output/"
    fi
    
    if [ -f "$outpath/NAVIP_First_Output/All_VCF.vcf" ] && [ -f "$outpath/NAVIP_Second_Output/All_VCF.vcf" ]; then
      if cmp -s "$outpath/NAVIP_First_Output/All_VCF.vcf" "$outpath/NAVIP_Second_Output/All_VCF.vcf"; then
        if rmkdir "$outpath/SFA_Output"; then
          rm -rf "$outpath/SFA_First_Output" "$outpath/SFA_Second_Output"
          python3 "$(dirname "$0")/navip.py" --mode sfa --innavipvcf "$outpath/NAVIP_First_Output/All_VCF.vcf" --outpath "$outpath/SFA_Output/"
        fi
      else
        if rmkdirs "$outpath/SFA_First_Output" "$outpath/SFA_Second_Output"; then
          rm -rf "$outpath/SFA_Output"
          python3 "$(dirname "$0")/navip.py" --mode sfa --innavipvcf "$outpath/NAVIP_First_Output/All_VCF.vcf" --outpath "$outpath/SFA_First_Output/"
          python3 "$(dirname "$0")/navip.py" --mode sfa --innavipvcf "$outpath/NAVIP_Second_Output/All_VCF.vcf" --outpath "$outpath/SFA_Second_Output/"
        fi
      fi
    fi
    
  fi
fi
