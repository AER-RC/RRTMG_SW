#!/bin/bash
exe=$(ls ../build/rrtmg_sw*)
platform=${exe#"../build/rrtmg_sw_v5.00_"}
outfld="./outputs_for_tests/${platform}"
allow_fail=$1
if [[ ${allow_fail} = "--allow-failed" ]] ; then
      af="allow"
      echo "Failed tests do not exit from testing $af"
fi
if [[ ! ${af} == "allow" ]] ; then
    echo "Testing will exit on the first failed test"
fi

declare -a test_input_mls_cld=("cld-imca0-icld2" "cld-imca1-icld2" "cld-imca1-icld2" \
                               "cld-imca1-icld5-idcor0" "cld-imca1-icld5-idcor1")


declare -a test_in_mls_cld=("cld_rrtm-cld6"   "cld_rrtm-cld5"    "cld_rrtm-cld7" \
                            "cld_rrtm-cld7" "cld_rrtm-cld7")

declare -a test_out_mls_cld=("cld6-imca0-icld2" "cld5-imca1-icld2" "cld7-imca1-icld2" \
                              "cld7-imca1-icld5-idcor0" "cld7-imca1-icld5-idcor1")

declare -a test_input_mls_clr=("clr" "clr-aer12" "clr-sza45-isolvar0_tsi_avg" \
                        "clr-sza45-isolvar1_tsi_max" "clr-sza45-isolvar1_tsi_min" \
                        "clr-sza45-isolvar2_01Jan1950" "clr-sza45-isolvar3_bndscl_tsi_max"
                        "clr-sza65" )
declare -a test_input_others=("MLW-clr" "SAW-clr" "TROP-clr")

    cd  ../test/

echo "testing MLS-cld:"

fail=""
for i in "${!test_input_mls_cld[@]}" ; do
    file="input_rrtm_MLS-${test_input_mls_cld[$i]}"
    file_in="in_${test_in_mls_cld[$i]}"
    file_out="${outfld}/output_rrtm_MLS-${test_out_mls_cld[$i]}"
    echo "testing ${file} with ${file_in}"
    cp -f "./inputs_for_tests/${file}" INPUT_RRTM
    cp -f "./inputs_for_tests/${file_in}" IN_CLD_RRTM
    ./rrtmg_sw
    if [[ -f OUTPUT_RRTM ]] ; then
        check=$(diff -q OUTPUT_RRTM ${file_out})
        if [[ ! "${check}" == "" ]] ; then
            echo "FAIL ${file} with ${file_in}"
            echo $(diff OUTPUT_RRTM ${file_out})
            if [[ ! ${af} == "allow" ]] ; then
              fail=1
              break
            fi
            rm -f INPUT_RRTM OUTPUT_RRTM IN_CLD_RRTM
        else
            echo "PASS ${file} with ${file_in}"
            rm -f INPUT_RRTM OUTPUT_RRTM IN_CLD_RRTM
        fi
    else
        echo "FAIL ${file} with ${file_in}"
        echo "rrtmg_sw did not execute properly"
        fail="1"
        break
    fi
done

if [[ $fail == "1" ]] ; then
    exit 1
fi


echo "testing MLS-clr:"
for i in "${!test_input_mls_clr[@]}" ; do
    file="input_rrtm_MLS-${test_input_mls_clr[$i]}"
    file_out="${outfld}/output_rrtm_MLS-${test_input_mls_clr[$i]}"
    echo "testing ${file}"
    cp -f "./inputs_for_tests/${file}" INPUT_RRTM
    if [[ ${test_input_mls_clr[$i]} == "clr-aer12" ]] ; then
        cp -f ./inputs_for_tests/in_aer_rrtm-aer12 IN_AER_RRTM
    fi
    ./rrtmg_sw
    diff -q OUTPUT_RRTM ${file_out}
    if [[ -f OUTPUT_RRTM ]] ; then
        check=$(diff -q OUTPUT_RRTM ${file_out})
        if [[ ! "${check}" == "" ]] ; then
            echo "FAIL ${file}"
            echo $(diff OUTPUT_RRTM ${file_out})
            if [[ ! ${af} == "allow" ]] ; then
              fail="1"
              break
            fi
            rm -f INPUT_RRTM OUTPUT_RRTM
        else
            echo "PASS ${file}"
            rm -f INPUT_RRTM OUTPUT_RRTM
        fi
    else
        echo "FAIL ${file}"
        echo "rrtmg_sw did not execute properly"
        fail="1"
        break
    fi
    if [[ ${test_input_mls_clr[$i]} == "clr-aer12" ]] ; then
        rm -f IN_AER_RRTM
    fi
done

if [[ $fail == "1" ]] ; then
    exit 1
fi
rm  -f tape6 TAPE7
echo "testing others:"

for i in "${!test_input_others[@]}" ; do
    file="input_rrtm_${test_input_others[$i]}"
    file_out="${outfld}/output_rrtm_${test_input_others[$i]}"
    echo "testing ${file}"
    cp -f "./inputs_for_tests/${file}" INPUT_RRTM
    ./rrtmg_sw
    diff -q OUTPUT_RRTM ${file_out}
    if [[ -f OUTPUT_RRTM ]] ; then
        check=$(diff -q OUTPUT_RRTM ${file_out})
        if [[ ! "${check}" == "" ]] ; then
            echo "FAIL ${file}"
            echo $(diff OUTPUT_RRTM ${file_out})
            if [[ ! ${af} == "allow" ]] ; then
              fail=1
              break
            fi
        else
            echo "PASS ${file}"
            rm -f INPUT_RRTM IN_CLD_RRTM OUTPUT_RRTM
        fi
    else
        echo "FAIL ${file}"
        echo "rrtmg_sw did not execute properly"
        fail="1"
        break
    fi
done

if [[ $fail == "1" ]] ; then
    exit 1
fi