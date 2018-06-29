
#ifndef _PARSE_INPUT_H_
#define _PARSE_INPUT_H_

#include <map>
#include <string>

namespace sgUtil {

   typedef std::map< std::string, std::string > ParseInputReturn_t;
   typedef std::map< std::string, std::string >::iterator ParseInputReturn_iterator_t;
   typedef std::map< std::string, std::string >::const_iterator ParseInputReturn_const_iterator_t;

   int ParseInput( const char* input_filename, ParseInputReturn_t& input_map );

};

#endif
