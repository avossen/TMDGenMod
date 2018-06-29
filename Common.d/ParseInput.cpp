
#include "ParseInput.h"

#include <map>

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>

using std::cerr;
using std::endl;

namespace sgUtil {

   int ParseInput( const char* input_filename, ParseInputReturn_t& input_map ){

      input_map.clear();

      if( !input_filename ){
         cerr << "ERROR: ParseInput given null pointer to file name" << endl;
         return 127;
      };

      std::ifstream fin( input_filename );
      if( !fin ){
         cerr << "ERROR opening file '" << input_filename << "'" << endl;
         return 127;
      };

      std::string line, key, obj, stemp;

      fin.peek();
      while( !fin.eof() ){
         obj = "";
         getline( fin, line );

         // trim the white space
         size_t pos = line.find_first_not_of(" \t");

         if( line.size() && !fin.eof() && pos != std::string::npos && line[pos] != '#' && line[pos] != '!'){
            std::istringstream ss(line);
            ss >> key;

            //cerr << "Found key " << key << endl;

            if( key.size() ){
               ss >> obj;
               ss >> stemp;

               while( ss ){
                  obj += ' ';
                  obj += stemp;
                  ss >> stemp;
               };

               //cerr << "\tobj = " << obj << endl;

               if( obj.size() ){
                  // remove '!' and anything following
                  size_t pos = obj.find_first_of('!');
                  if( pos != std::string::npos ){
                     size_t pos2 = obj.find_last_not_of(" \t", pos-1, 2 );

                     if( pos2 == std::string::npos ){
                        obj = obj.substr( 0, pos );
                     } else {
                        obj = obj.substr( 0, pos2+1 );
                     };
                  };

                  //cerr << "\tobj = " << obj << endl;

                  std::pair< ParseInputReturn_iterator_t, bool > insert_result = input_map.insert( std::pair< std::string, std::string >( key, obj ) );
                  if( !insert_result.second ){
                     // something already exists with this key, so delete it and add this one
                     input_map.erase( insert_result.first );

                     // can't fail this time, as it can have only one of each key type, and it was just removed
                     input_map.insert( std::pair< std::string, std::string >( key, obj ) );
                  };
               };
            };
         };

         fin.peek();
      };

//       ParseInputReturn_iterator_t it;
//       for( it = input_map.begin(); it != input_map.end(); ++it )
//          cerr << "***\t" << it->first << ' ' << it->second << endl;


      return 0;
   };

};
