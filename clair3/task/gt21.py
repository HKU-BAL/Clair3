from enum import IntEnum

GT21_LABELS = [
    'AA',
    'AC',
    'AG',
    'AT',
    'CC',
    'CG',
    'CT',
    'GG',
    'GT',
    'TT',
    'DelDel',
    'ADel',
    'CDel',
    'GDel',
    'TDel',
    'InsIns',
    'AIns',
    'CIns',
    'GIns',
    'TIns',
    'InsDel'
]
GT21_LABELS_MAP = dict(zip(GT21_LABELS, range(0, 21)))


class GT21_Type(IntEnum):
    AA = 0
    AC = 1
    AG = 2
    AT = 3
    CC = 4
    CG = 5
    CT = 6
    GG = 7
    GT = 8
    TT = 9
    DelDel = 10
    ADel = 11
    CDel = 12
    GDel = 13
    TDel = 14
    InsIns = 15
    AIns = 16
    CIns = 17
    GIns = 18
    TIns = 19
    InsDel = 20


def gt21_label_from(gt21_enum):
    try:
        return GT21_LABELS[gt21_enum]
    except:
        return ""


def gt21_enum_from_label(gt21_label):
    return GT21_LABELS_MAP[gt21_label]


def partial_label_from(ref, alt):
    if len(ref) > len(alt):
        return "Del"
    elif len(ref) < len(alt):
        return "Ins"
    return alt[0]


def mix_two_partial_labels(label1, label2):
    # AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
    if len(label1) == 1 and len(label2) == 1:
        return label1 + label2 if label1 <= label2 else label2 + label1

    # ADel, CDel, GDel, TDel, AIns, CIns, GIns, TIns
    tmp_label1, tmp_label2 = label1, label2
    if len(label1) > 1 and len(label2) == 1:
        tmp_label1, tmp_label2 = label2, label1
    if len(tmp_label2) > 1 and len(tmp_label1) == 1:
        return tmp_label1 + tmp_label2

    # InsIns, DelDel
    if len(label1) > 0 and len(label2) > 0 and label1 == label2:
        return label1 + label2

    # InsDel
    return gt21_label_from(GT21_Type.InsDel)


def gt21_enum_from(reference, alternate, genotype_1, genotype_2, alternate_arr=None):
    if alternate_arr is not None:
        partial_labels = [partial_label_from(reference, alternate) for alternate in alternate_arr]
        gt21_label = mix_two_partial_labels(partial_labels[0], partial_labels[1])
        return gt21_enum_from_label(gt21_label)

    alternate_arr = alternate.split(',')
    if len(alternate_arr) == 1:
        alternate_arr = (
            [reference if genotype_1 == 0 or genotype_2 == 0 else alternate_arr[0]] +
            alternate_arr
        )

    partial_labels = [partial_label_from(reference, alternate) for alternate in alternate_arr]
    gt21_label = mix_two_partial_labels(partial_labels[0], partial_labels[1])

    return gt21_enum_from_label(gt21_label)


HOMO_SNP_GT21 = [GT21_Type.AA, GT21_Type.CC, GT21_Type.GG, GT21_Type.TT]
HOMO_SNP_LABELS = [gt21_label_from(gt21_enum) for gt21_enum in HOMO_SNP_GT21]

HETERO_SNP_GT21 = [GT21_Type.AC, GT21_Type.AG, GT21_Type.AT, GT21_Type.CG, GT21_Type.CT, GT21_Type.GT]
HETERO_SNP_LABELS = [gt21_label_from(gt21_enum) for gt21_enum in HETERO_SNP_GT21]
