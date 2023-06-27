       for i in {1..22} ; do \
                V2=$(bc -l <<< "0.1*e((${i}-1)*l(1.5))")
                echo "LINKING JOB NR.${i} FOR V2=${V2}"
        done

