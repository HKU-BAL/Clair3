import sys
from importlib import import_module
from shared.param_p import REPO_NAME

DATA_PREP_SCRIPTS_FOLDER="preprocess"
DEEP_LEARNING_FOLDER="clair3"
POST_PROCESS_SCRIPTS_FOLDER="clair3.metrics"

deep_learning_folder = [
    "CallVarBam",
    "CallVariants",
    "Train",
]

data_preprocess_folder = [
    "GetTruth",
    "Tensor2Bin",
    'RealignReads',
    'CreateTensorPileup',
    "CreateTensorFullAlignment",
    'CreateTrainingTensor',
    'SplitExtendBed',
    'MergeBin',
    'MergeVcf',
    'SelectHetSnp',
    'SelectCandidates',
    'UnifyRepresentation',
    'CheckEnvs',
    'SortVcf',
    'SelectQual'
]

post_process_scripts_folder = [
    'GetOverallMetrics',
]

def directory_for(submodule_name):
    if submodule_name in deep_learning_folder:
        return DEEP_LEARNING_FOLDER
    if submodule_name in data_preprocess_folder:
        return DATA_PREP_SCRIPTS_FOLDER
    if submodule_name in post_process_scripts_folder:
        return POST_PROCESS_SCRIPTS_FOLDER
    return ""


def print_help_messages():
    from textwrap import dedent
    print(dedent("""\
        {0} submodule invocator:
            Usage: python clair3.py [submodule] [options of the submodule]
        Available data preparation submodules:\n{1}
        Available clair submodules:\n{2}
        Available post processing submodules:\n{3}
        """.format(
            REPO_NAME,
            "\n".join("          - %s" % submodule_name for submodule_name in data_preprocess_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in deep_learning_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in post_process_scripts_folder),
        )
    ))


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print_help_messages()
        sys.exit(0)

    submodule_name = sys.argv[1]
    if (
        submodule_name not in deep_learning_folder and
        submodule_name not in data_preprocess_folder and
        submodule_name not in post_process_scripts_folder
    ):
        sys.exit("[ERROR] Submodule %s not found." % (submodule_name))

    directory = directory_for(submodule_name)
    submodule = import_module("%s.%s" % (directory, submodule_name))

    # filter arguments (i.e. filter clair3.py) and add ".py" for that submodule
    sys.argv = sys.argv[1:]
    sys.argv[0] += (".py")

    # Note: need to make sure every submodule contains main() method
    submodule.main()

    sys.exit(0)


if __name__ == "__main__":
    main()
