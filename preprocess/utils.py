

import mpmath as math
from cffi import FFI
import logging


class mathcalculator(object):


    def __init__(self,speedUp=True):

        self.LOG_10 = math.log(10.0)
        self.maxPhredScore = 255
        self.speedUp = speedUp
        if(speedUp):
            self._creatCFFIFunc()
            
            '''
            print("speed up!")
            self.ffi = FFI()
            self.ffi.cdef("""double log10p_to_phred(double log10p);""")
            self.lib = self.ffi.verify(sources=['preprocess/mathCalculator/clairMath.c'])
            '''

    
    def _creatCFFIFunc(self):
        self.ffi = FFI()
        self.ffi.cdef("""
                double log10p_to_phred(double log10p);
                double log10sumexp(double log10_array[],int n_array);
                double getMyMaxItem(double list[],int n_list);
                """)
        self.lib = self.ffi.verify("""
                        #include<math.h>
                        #include<stdio.h>
                        double LOG_10 = log(10.0);
                        double log10p_to_phred(double log10p){
                            double ptrue;
                            ptrue = pow(10,log10p);
                            if(ptrue==1){
                                return 255;
                            }

                            return -10*(log(1.0-ptrue)/LOG_10);
                        }
                       double f_log10(double myInput){
                            double res;

                            res=log(myInput)/LOG_10;
                            return res;
                       }
                       double getMyMaxItem(double list[],int n_list){
                           double curMax;
                           int i;

                           curMax = list[0];
                           for(i=1;i<=n_list;i++){
                               if(list[i]>curMax){
                                   curMax= list[i];
                               }
                            }


                           return curMax;
                       } 


                       double log10sumexp(double log10_array[],int n_array){

                           double m, mySum,tmp;
                           int i;
        
                           m = getMyMaxItem(log10_array,n_array);
                           mySum = 0.0;
                           for(i=0;i<n_array;i++){

                               tmp = pow(10,log10_array[i]-m);
                               mySum += tmp;
                           }

                           return log(mySum)/LOG_10;

                       }



                       """
                       )

          
    def log10p_to_phred(self,log10p):

        '''
        log10(p) to -10log10(1-p)

        '''
        #print('log10p',log10p)
        if(self.speedUp):
            ''' 
            res1 = self.lib.log10p_to_phred(log10p)
            ptrue = math.power(10,log10p)
            if(ptrue==1):
                res2 = 255
            else:

                res2 = -10*(math.log1p(-ptrue)/self.LOG_10)
            logging.info('log10p_to_phred c version '+str(res1))
            logging.info('log10p_to_phred pypy version '+str(res2))
            '''
            return round(self.lib.log10p_to_phred(log10p),6)
        
        ptrue = math.power(10,log10p)
        
        if(ptrue==1):
            return 255
        #print('ptrue',ptrue,'math.log1p(ptrue)/self.LOG_10',math.log1p(-ptrue)/self.LOG_10,'return value',-10*(math.log1p(-ptrue)/self.LOG_10))
        return round(-10*(math.log1p(-ptrue)/self.LOG_10),6)
        
        
    
    def normalizeLog10Likelihood(self,log10prob):
        
        

        pass
    def log10sumexp(self,log10_array):

        '''

        return value approxiameting to log10(sum(10^log10_probs))
        '''

        if(self.speedUp):
            
            n_array = self.ffi.cast("int", len(log10_array))
            log10_array_c = self.ffi.new("double []",log10_array)
            '''
            res1 = self.lib.log10sumexp(log10_array_c,n_array)
            m = max(log10_array)
            res2 = m + math.log10(sum(pow(10.0, x - m) for x in log10_array))
            logging.info("log10sumexp c version "+str(res1))
            logging.info("log10sumexp pypy version "+str(res2))
            '''
            return self.lib.log10sumexp(log10_array_c,n_array) 
            
        m = max(log10_array)
        return  m + math.log10(sum(pow(10.0, x - m) for x in log10_array))
   

    def normalize_log10_prob(self,log10_prob):

        '''
        normalize the log10 probs that make sum of the probabilities close to 1
    
        '''
        # keep 6 digits
        lse = round(self.log10sumexp(log10_prob),6)
        #print('lse',lse)
        normalized_log10_probs = [] 
        for x in log10_prob:
            normalized_log10_probs.append(min(x-lse,0))
   
        return normalized_log10_probs











