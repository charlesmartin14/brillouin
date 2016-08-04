#!/usr/bin/env  bash
cat $1  | sed -e 's/ ∗ ∗/**/g' | sed -e 's/∗ ∗/**/g' | tr '∗' '*' | tr '−' '-' | sed -e 's/ d0/d0/g' | sed -e 's/ (/(/g' | sed -e 's/ )/)/g' | sed -e 's/( /(/g' | sed -e 's/ )/)/g' | sed -e 's/ ,/,/g' | sed -e 's/, /,/g'  | sed -e 's/ \. /./g' | sed -e 's/^ [0-9][0-9] /    /' | sed -e 's/^ [0-9][0-9][0-9]/    /' 


