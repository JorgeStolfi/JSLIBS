#! /bin/sed -f
# Last edited on 2006-04-15 11:15:47 by stolfi

# Converts clients of "SPParams.h" into clients of "argparser.h"

s/SPParams_Error/argparser_error/g
s/SPParams_Finish/argparser_finish/g
s/SPParams_GetIntList/argparser_get_int_list/g
s/SPParams_GetKeyword/argparser_get_keyword/g
s/SPParams_GetNextDir/argparser_get_next_dir/g
s/SPParams_GetNextDouble/argparser_get_next_double/g
s/SPParams_GetNextInt/argparser_get_next_int/g
s/SPParams_GetNextR3/argparser_get_next_r3/g
s/SPParams_GetNextR4/argparser_get_next_r4/g
s/SPParams_GetNext/argparser_get_next/g
s/SPParams_IsNext/argparser_is_next/g
s/SPParams_KeywordPresent/argparser_keyword_present/g
s/SPParams_MatchNext/argparser_get_next_keyword/g
s/SPParams_NewT/argparser_new/g
s/SPParams_SetUsage/argparser_set_help/g
s/SPParams_SkipParsed/argparser_skip_parsed/g
s/SPParams_TestNext/argparser_keyword_present_next/g
s/SPParams_T/argparser_t/g
s/SPParams/argparser/g
