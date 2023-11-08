#ifndef _GEN_PAR_CPP_H_
#define _GEN_PAR_CPP_H_

const int CHECK_DEGREE[192]=
{
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

const int LINK_NODES[768]=
{
   72, 144, 216, 360, /* [ 0] */
   73, 145, 217, 361, /* [ 1] */
   74, 146, 218, 362, /* [ 2] */
   75, 147, 219, 363, /* [ 3] */
   76, 148, 220, 364, /* [ 4] */
   77, 149, 221, 365, /* [ 5] */
   78, 150, 222, 366, /* [ 6] */
   79, 151, 223, 367, /* [ 7] */
   80, 152, 224, 368, /* [ 8] */
   81, 153, 225, 369, /* [ 9] */
   82, 154, 226, 370, /* [10] */
   83, 155, 227, 371, /* [11] */
   84, 156, 228, 372, /* [12] */
   85, 157, 229, 373, /* [13] */
   86, 158, 230, 374, /* [14] */
   87, 159, 231, 375, /* [15] */
   88, 160, 232, 376, /* [16] */
   89, 161, 233, 377, /* [17] */
   90, 162, 234, 378, /* [18] */
   91, 163, 235, 379, /* [19] */
   92, 164, 236, 380, /* [20] */
   93, 165, 237, 381, /* [21] */
   94, 166, 238, 382, /* [22] */
   95, 167, 239, 383, /* [23] */
   72,  96, 192, 336, /* [24] */
   73,  97, 193, 337, /* [25] */
   74,  98, 194, 338, /* [26] */
   75,  99, 195, 339, /* [27] */
   76, 100, 196, 340, /* [28] */
   77, 101, 197, 341, /* [29] */
   78, 102, 198, 342, /* [30] */
   79, 103, 199, 343, /* [31] */
   80, 104, 200, 344, /* [32] */
   81, 105, 201, 345, /* [33] */
   82, 106, 202, 346, /* [34] */
   83, 107, 203, 347, /* [35] */
   84, 108, 204, 348, /* [36] */
   85, 109, 205, 349, /* [37] */
   86, 110, 206, 350, /* [38] */
   87, 111, 207, 351, /* [39] */
   88, 112, 208, 352, /* [40] */
   89, 113, 209, 353, /* [41] */
   90, 114, 210, 354, /* [42] */
   91, 115, 211, 355, /* [43] */
   92, 116, 212, 356, /* [44] */
   93, 117, 213, 357, /* [45] */
   94, 118, 214, 358, /* [46] */
   95, 119, 215, 359, /* [47] */
   24, 120, 152, 312, /* [48] */
   25, 121, 153, 313, /* [49] */
   26, 122, 154, 314, /* [50] */
   27, 123, 155, 315, /* [51] */
   28, 124, 156, 316, /* [52] */
   29, 125, 157, 317, /* [53] */
   30, 126, 158, 318, /* [54] */
   31, 127, 159, 319, /* [55] */
   32, 128, 160, 320, /* [56] */
   33, 129, 161, 321, /* [57] */
   34, 130, 162, 322, /* [58] */
   35, 131, 163, 323, /* [59] */
   36, 132, 164, 324, /* [60] */
   37, 133, 165, 325, /* [61] */
   38, 134, 166, 326, /* [62] */
   39, 135, 167, 327, /* [63] */
   40, 136, 144, 328, /* [64] */
   41, 137, 145, 329, /* [65] */
   42, 138, 146, 330, /* [66] */
   43, 139, 147, 331, /* [67] */
   44, 140, 148, 332, /* [68] */
   45, 141, 149, 333, /* [69] */
   46, 142, 150, 334, /* [70] */
   47, 143, 151, 335, /* [71] */
    0, 168, 240, 362, /* [72] */
    1, 169, 241, 363, /* [73] */
    2, 170, 242, 364, /* [74] */
    3, 171, 243, 365, /* [75] */
    4, 172, 244, 366, /* [76] */
    5, 173, 245, 367, /* [77] */
    6, 174, 246, 368, /* [78] */
    7, 175, 247, 369, /* [79] */
    8, 176, 248, 370, /* [80] */
    9, 177, 249, 371, /* [81] */
   10, 178, 250, 372, /* [82] */
   11, 179, 251, 373, /* [83] */
   12, 180, 252, 374, /* [84] */
   13, 181, 253, 375, /* [85] */
   14, 182, 254, 376, /* [86] */
   15, 183, 255, 377, /* [87] */
   16, 184, 256, 378, /* [88] */
   17, 185, 257, 379, /* [89] */
   18, 186, 258, 380, /* [90] */
   19, 187, 259, 381, /* [91] */
   20, 188, 260, 382, /* [92] */
   21, 189, 261, 383, /* [93] */
   22, 190, 262, 360, /* [94] */
   23, 191, 263, 361, /* [95] */
   48, 205, 245, 315, /* [96] */
   49, 206, 246, 316, /* [97] */
   50, 207, 247, 317, /* [98] */
   51, 208, 248, 318, /* [99] */
   52, 209, 249, 319, /* [100] */
   53, 210, 250, 320, /* [101] */
   54, 211, 251, 321, /* [102] */
   55, 212, 252, 322, /* [103] */
   56, 213, 253, 323, /* [104] */
   57, 214, 254, 324, /* [105] */
   58, 215, 255, 325, /* [106] */
   59, 192, 256, 326, /* [107] */
   60, 193, 257, 327, /* [108] */
   61, 194, 258, 328, /* [109] */
   62, 195, 259, 329, /* [110] */
   63, 196, 260, 330, /* [111] */
   64, 197, 261, 331, /* [112] */
   65, 198, 262, 332, /* [113] */
   66, 199, 263, 333, /* [114] */
   67, 200, 240, 334, /* [115] */
   68, 201, 241, 335, /* [116] */
   69, 202, 242, 312, /* [117] */
   70, 203, 243, 313, /* [118] */
   71, 204, 244, 314, /* [119] */
    0, 133, 264, 348, /* [120] */
    1, 134, 265, 349, /* [121] */
    2, 135, 266, 350, /* [122] */
    3, 136, 267, 351, /* [123] */
    4, 137, 268, 352, /* [124] */
    5, 138, 269, 353, /* [125] */
    6, 139, 270, 354, /* [126] */
    7, 140, 271, 355, /* [127] */
    8, 141, 272, 356, /* [128] */
    9, 142, 273, 357, /* [129] */
   10, 143, 274, 358, /* [130] */
   11, 120, 275, 359, /* [131] */
   12, 121, 276, 336, /* [132] */
   13, 122, 277, 337, /* [133] */
   14, 123, 278, 338, /* [134] */
   15, 124, 279, 339, /* [135] */
   16, 125, 280, 340, /* [136] */
   17, 126, 281, 341, /* [137] */
   18, 127, 282, 342, /* [138] */
   19, 128, 283, 343, /* [139] */
   20, 129, 284, 344, /* [140] */
   21, 130, 285, 345, /* [141] */
   22, 131, 286, 346, /* [142] */
   23, 132, 287, 347, /* [143] */
   24,  97, 180, 288, /* [144] */
   25,  98, 181, 289, /* [145] */
   26,  99, 182, 290, /* [146] */
   27, 100, 183, 291, /* [147] */
   28, 101, 184, 292, /* [148] */
   29, 102, 185, 293, /* [149] */
   30, 103, 186, 294, /* [150] */
   31, 104, 187, 295, /* [151] */
   32, 105, 188, 296, /* [152] */
   33, 106, 189, 297, /* [153] */
   34, 107, 190, 298, /* [154] */
   35, 108, 191, 299, /* [155] */
   36, 109, 168, 300, /* [156] */
   37, 110, 169, 301, /* [157] */
   38, 111, 170, 302, /* [158] */
   39, 112, 171, 303, /* [159] */
   40, 113, 172, 304, /* [160] */
   41, 114, 173, 305, /* [161] */
   42, 115, 174, 306, /* [162] */
   43, 116, 175, 307, /* [163] */
   44, 117, 176, 308, /* [164] */
   45, 118, 177, 309, /* [165] */
   46, 119, 178, 310, /* [166] */
   47,  96, 179, 311, /* [167] */
   48, 226, 282, 296, /* [168] */
   49, 227, 283, 297, /* [169] */
   50, 228, 284, 298, /* [170] */
   51, 229, 285, 299, /* [171] */
   52, 230, 286, 300, /* [172] */
   53, 231, 287, 301, /* [173] */
   54, 232, 264, 302, /* [174] */
   55, 233, 265, 303, /* [175] */
   56, 234, 266, 304, /* [176] */
   57, 235, 267, 305, /* [177] */
   58, 236, 268, 306, /* [178] */
   59, 237, 269, 307, /* [179] */
   60, 238, 270, 308, /* [180] */
   61, 239, 271, 309, /* [181] */
   62, 216, 272, 310, /* [182] */
   63, 217, 273, 311, /* [183] */
   64, 218, 274, 288, /* [184] */
   65, 219, 275, 289, /* [185] */
   66, 220, 276, 290, /* [186] */
   67, 221, 277, 291, /* [187] */
   68, 222, 278, 292, /* [188] */
   69, 223, 279, 293, /* [189] */
   70, 224, 280, 294, /* [190] */
   71, 225, 281, 295  /* [191] */
};
const int ALPHA[768]=
{
   60,  32,  23,  45, /* [ 0] */
   33,   5,  18,  59, /* [ 1] */
   16,  29,   7,  44, /* [ 2] */
   10,  19,  47,  32, /* [ 3] */
   22,  37,   9,  63, /* [ 4] */
   27,  40,  18,  55, /* [ 5] */
   41,  63,  50,  15, /* [ 6] */
   47,  12,  38,  60, /* [ 7] */
   41,  19,  28,  56, /* [ 8] */
   28,  50,   2,  37, /* [ 9] */
   34,  47,  62,  25, /* [10] */
   29,  44,  16,   7, /* [11] */
   36,   8,  21,  62, /* [12] */
   42,  55,   7,  33, /* [13] */
   11,  39,   2,  24, /* [14] */
   34,  12,  21,  49, /* [15] */
    9,  24,  50,  59, /* [16] */
   13,  22,  50,  35, /* [17] */
   40,  27,  55,  18, /* [18] */
   15,  56,  30,   2, /* [19] */
   46,  55,   5,  20, /* [20] */
   15,  52,  24,  37, /* [21] */
   62,  25,  34,  47, /* [22] */
   18,  46,  31,   9, /* [23] */
   37,  22,  63,   9, /* [24] */
   50,  35,  22,  13, /* [25] */
   16,  29,   7,  44, /* [26] */
   38,  51,  29,   3, /* [27] */
    7,  48,  57,  22, /* [28] */
   55,   1,  29,  14, /* [29] */
   17,   8,  45,  30, /* [30] */
   48,   7,  22,  57, /* [31] */
   21,  47,   6,  56, /* [32] */
   24,  59,  50,   9, /* [33] */
   39,  52,  30,   4, /* [34] */
    6,  21,  47,  56, /* [35] */
   48,  61,  39,  13, /* [36] */
   51,  23,  14,  36, /* [37] */
   55,  14,  29,   1, /* [38] */
   34,   6,  19,  60, /* [39] */
   60,  12,  38,  47, /* [40] */
   56,   6,  21,  47, /* [41] */
   31,  18,  46,   9, /* [42] */
   29,  14,  55,   1, /* [43] */
    3,  53,  18,  44, /* [44] */
   33,  24,  61,  46, /* [45] */
   33,   7,  42,  55, /* [46] */
    6,  60,  19,  34, /* [47] */
   58,  23,  49,   8, /* [48] */
   27,  62,  12,  53, /* [49] */
   32,  23,  60,  45, /* [50] */
   63,  22,   9,  37, /* [51] */
   61,   7,  20,  35, /* [52] */
   56,  47,  21,   6, /* [53] */
   39,  26,  17,  54, /* [54] */
   24,  37,  15,  52, /* [55] */
   14,  36,  23,  51, /* [56] */
   27,  36,   1,  49, /* [57] */
   62,  25,  47,  34, /* [58] */
   27,   5,  42,  14, /* [59] */
   59,   9,  50,  24, /* [60] */
   48,  33,  11,  20, /* [61] */
   14,  40,  62,  49, /* [62] */
   45,  10,  58,  36, /* [63] */
   50,   2,  37,  28, /* [64] */
   13,  28,  54,  63, /* [65] */
   48,  39,  13,  61, /* [66] */
   47,  38,  60,  12, /* [67] */
   55,  27,  40,  18, /* [68] */
   42,   5,  14,  27, /* [69] */
   58,  32,  17,   4, /* [70] */
   38,  53,  16,  25, /* [71] */
    1,  55,  29,  14, /* [72] */
   50,  37,   2,  28, /* [73] */
   14,   5,  27,  42, /* [74] */
   27,  14,  42,   5, /* [75] */
   17,  26,  54,  39, /* [76] */
   12,  53,  27,  62, /* [77] */
   27,  40,  18,  55, /* [78] */
   22,  50,  35,  13, /* [79] */
   31,  22,  44,  59, /* [80] */
   29,   3,  51,  38, /* [81] */
   19,  45,   4,  54, /* [82] */
   10,  51,  25,  60, /* [83] */
   48,  13,  61,  39, /* [84] */
   29,   1,  55,  14, /* [85] */
   17,  30,  45,   8, /* [86] */
   27,  18,  55,  40, /* [87] */
   59,  18,  33,   5, /* [88] */
   46,  24,  61,  33, /* [89] */
   63,  54,  28,  13, /* [90] */
   12,  47,  38,  60, /* [91] */
   23,  32,  45,  60, /* [92] */
   12,  21,  34,  49, /* [93] */
    6,  15,  43,  28, /* [94] */
   15,  52,  24,  37, /* [95] */
   39,  61,  13,  48, /* [96] */
   20,  46,   5,  55, /* [97] */
   32,  17,  58,   4, /* [98] */
   51,  14,  23,  36, /* [99] */
   14,  27,   5,  42, /* [100] */
   47,  60,  12,  38, /* [101] */
   24,  15,  52,  37, /* [102] */
   30,   2,  56,  15, /* [103] */
   11,  61,  26,  52, /* [104] */
   62,   8,  21,  36, /* [105] */
   13,  54,  28,  63, /* [106] */
    7,  48,  22,  57, /* [107] */
   19,  34,   6,  60, /* [108] */
   33,  42,  55,   7, /* [109] */
   51,  25,  10,  60, /* [110] */
    4,  45,  54,  19, /* [111] */
   54,  13,  63,  28, /* [112] */
   45,  54,   4,  19, /* [113] */
   22,  48,   7,  57, /* [114] */
   30,  17,   8,  45, /* [115] */
   45,  32,  23,  60, /* [116] */
   51,  38,   3,  29, /* [117] */
    6,  43,  28,  15, /* [118] */
   33,  11,  48,  20, /* [119] */
   44,  59,  31,  22, /* [120] */
   57,  44,   9,  35, /* [121] */
   38,  25,  53,  16, /* [122] */
    6,  47,  21,  56, /* [123] */
   23,  10,   1,  38, /* [124] */
   56,   6,  21,  47, /* [125] */
   35,   9,  57,  44, /* [126] */
   26,   4,  41,  13, /* [127] */
   62,  53,  27,  12, /* [128] */
   39,  61,  13,  48, /* [129] */
   50,  13,  22,  35, /* [130] */
    2,  39,  24,  11, /* [131] */
   21,  12,  49,  34, /* [132] */
   19,  56,  41,  28, /* [133] */
   51,  16,   1,  42, /* [134] */
   30,   8,  17,  45, /* [135] */
   27,  42,  14,   5, /* [136] */
    8,  34,  43,  56, /* [137] */
   17,  39,  26,  54, /* [138] */
   30,  56,   2,  15, /* [139] */
   14,  62,  49,  40, /* [140] */
   30,  43,  58,  21, /* [141] */
   14,   1,  55,  29, /* [142] */
   56,  19,  41,  28, /* [143] */
   17,  26,  54,  39, /* [144] */
   18,  53,  44,   3, /* [145] */
    1,  51,  42,  16, /* [146] */
   57,  16,   3,  31, /* [147] */
   27,  62,  12,  53, /* [148] */
   54,  17,  26,  39, /* [149] */
   40,  62,  49,  14, /* [150] */
   53,  16,  25,  38, /* [151] */
   51,  10,  60,  25, /* [152] */
   30,   2,  15,  56, /* [153] */
   12,  62,  53,  27, /* [154] */
   50,  59,  24,   9, /* [155] */
    7,  42,  55,  33, /* [156] */
   21,  43,  30,  58, /* [157] */
   35,  57,   9,  44, /* [158] */
   39,  17,  26,  54, /* [159] */
   62,  25,  47,  34, /* [160] */
   32,  23,  60,  45, /* [161] */
   40,  25,  12,   3, /* [162] */
   58,  49,   8,  23, /* [163] */
   22,   9,  37,  63, /* [164] */
   39,  61,  48,  13, /* [165] */
   47,  21,   6,  56, /* [166] */
   19,  47,  32,  10, /* [167] */
   33,  20,  11,  48, /* [168] */
   31,  22,  44,  59, /* [169] */
   51,  60,  10,  25, /* [170] */
    6,  56,  21,  47, /* [171] */
   37,  59,  46,  11, /* [172] */
   55,  27,  40,  18, /* [173] */
   24,   9,  59,  50, /* [174] */
    8,  56,  43,  34, /* [175] */
   44,  22,  59,  31, /* [176] */
   10,  45,  36,  58, /* [177] */
   31,   5,  40,  53, /* [178] */
   19,  41,  56,  28, /* [179] */
   15,  28,  43,   6, /* [180] */
   61,  11,  26,  52, /* [181] */
   57,   3,  16,  31, /* [182] */
   14,  27,   5,  42, /* [183] */
    5,  14,  27,  42, /* [184] */
   44,   7,  16,  29, /* [185] */
   49,  23,   8,  58, /* [186] */
   36,   1,  27,  49, /* [187] */
    9,  18,  46,  31, /* [188] */
   63,  37,  22,   9, /* [189] */
   43,  52,  17,   2, /* [190] */
    5,  40,  31,  53  /* [191] */
};

#endif
