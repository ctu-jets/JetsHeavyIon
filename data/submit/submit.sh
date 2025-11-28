#!/bin/bash
# -- baseFolder
baseFolder="$PWD"
# -- listOfFiles
listOfFiles="$PWD/$1"
# strip BASENAME from the path without the extension to get filelist_name range
filelist_name=$(basename "${listOfFiles}" .list)
# -- tree name
treeName="jetTree"
# -- root macro
rootMacro="runPicoHFJetMaker.C"
# -- bad run list file
badRunListFileName="BadRunList_14.list"
# -- production Id
productionId=$(date +%F)
# -- set STAR software version
starVersion="pro"
# --max number of files
maxNFiles=75

run_mode="$2"  # pass "embedding" or "data" from run.sh
isEmbedding="true"
if [[ "$run_mode" == "data" ]]; then
  isEmbedding="false"
fi


#================================================================
# -- submission xml file
templateXml="template.xml"

jobFolder="${baseFolder}/submit/${productionId}/job_${filelist_name}"
# -- job submission directory
mkdir -p "${jobFolder}"
cd "${jobFolder}"
# -- prepare folder
mkdir -p report err log list csh production

check=(
    ".sl73_x8664_gcc485" # dirs
    "StRoot"
    "StRoot/macros/${rootMacro}" # run macro
    "submit/${templateXml}"      # xml template
    "${badRunListFileName}"      # bad-run list
)

printf "Checking project … "
for item in "${check[@]}"; do
    path="${baseFolder}/${item}"
    [[ -e $path ]] || {
        echo "$item missing"
        exit 1
    }
    [[ -d $path || $item == *.xml || $item == "$badRunListFileName" ]] && ln -sf "$path"
done
[[ -e $listOfFiles ]] || {
    echo "$listOfFiles missing"
    exit 1
}
echo "ok"

# -- submit
generatedXml="generated_${filelist_name}.xml"

cat <<EOF >"${generatedXml}"
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE note [
<!ENTITY rootMacro "${rootMacro}">
<!ENTITY treeName "${treeName}">
<!ENTITY baseFolder "${baseFolder}">
<!ENTITY jobFolder "${jobFolder}">
<!ENTITY badRunListFileName "${badRunListFileName}">
<!ENTITY listOfFiles "${listOfFiles}">
<!ENTITY starVersion "${starVersion}">
<!ENTITY maxNFiles "${maxNFiles}">
<!ENTITY filelist_name "${filelist_name}">
<!ENTITY isEmbedding "${isEmbedding}">

]>

EOF
# -- add the rest of the xml file except the first line <?xml version="1.0" encoding="utf-8" ?>
tail -n +2 "${templateXml}" >>"${generatedXml}"

star-submit -u ie "${generatedXml}"
# star-submit-beta "${generatedXml}"
