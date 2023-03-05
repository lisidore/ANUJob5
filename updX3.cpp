#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]

//using namespace arma; I won't use this as sometimes function names can overlap



// [[Rcpp::export]]




Rcpp::List updX3(Rcpp::DataFrame events,int nrow_events,double edget,Rcpp::NumericVector dendistrt,Rcpp::NumericVector dendistrto,Rcpp::NumericVector nichedistrt,
                 Rcpp::NumericVector backdistrt, double a, double b,double u,int x,double sd, int N, double d) {  //not sure yet what this function will return
  
  
  
  int nrow_ev = nrow_events;
  // int total = sum(dendistrt);  ****unused, will be deleted
  Rcpp::NumericVector dendistrt2; //all of these variables need to be initialised
  Rcpp::NumericVector tmp1(dendistrto.size()); //initialise the vector
  Rcpp::IntegerVector temp_vec5;
  int accum = 0;
  
  Rcpp::IntegerVector events_type = events["type"];
  Rcpp::NumericVector events_n;
  Rcpp::NumericVector events_a = events["a"];
  Rcpp::NumericVector events_b = events["b"];
  Rcpp::IntegerVector events_xmin = events["xmin"];
  Rcpp::NumericVector events_time = events["time"];
  
  Rcpp::IntegerVector dendistrt_N(N); // will be used in the lines corresponding to R code: dendistrt <- numeric(N)
  
  //Start of declarations concerning R code ******** dendistrto <- dendistrto + c(rep(1,N),nichedistrt) ...
  int vect_size = N + nichedistrt.size();
  Rcpp::NumericVector rep_vec1 = Rcpp::rep(1.0,N);
  Rcpp::NumericVector rep_vec2;
  int w4 = 0;
  //End of declarations concerning R code ******** dendistrto <- dendistrto + c(rep(1,N),nichedistrt) ...
  //start of for loop to copy rep_vec1 into rep_vec2 in case the declaration and initialisation of rep_vec1 has limited it to a certain length
  
  for(int y = 0; y < N; y++){ 
    rep_vec2[y] = rep_vec1[y];
  } 
  
  //end of for loop to copy rep_vec1 into rep_vec2 
  
  double pt;
  //int tmp_abs = 0; //unused, will be deleted
  
  int temp1 = 0;
  int temp2 = 0;
  int temp3 = 0;
  int temp4 = 0;
  // Rcpp::IntegerVector mat_index;
  //  int vec_index = 0; //unused, will be deleted
  
  
  
  if (nrow_ev!=0 && sum(dendistrt)>0) { //start of top_level_if_1
    
    Rcpp::IntegerVector tmpo;
    
    if(events_time[0] > 0){ // start of top_level_if_2
      x = events_xmin[0];
      
      
      //****************************************************************************************************************   
      //equivalent to R code **** min(which(dendistrt>0))
      int tempo_index = -1;
      for(int v = 0; v < dendistrt.size(); v++){ //start of for 9
        if(dendistrt[v] > 0 && tempo_index < 0){
          tempo_index = v;  // minimum index for which dendistrt > 0 will be stored here
        } 
      } //end of for 9
      
      if(x > tempo_index) { //start of top_level_if_3
                dendistrt = dendistrt_N;      //correspoding to R code: dendistrt <- numeric(N) 
      }else{ //end of top_level_if_3, start of top_level_else_3
        
        if(x > 1){ //start of top_level_if_4
          
          //start of section corresponding to R code ******** if (x>1) {backdistrt[1:(x-1+N)] <- 0} 
          int vect_index1 = x - 1 + N;
          for(int q = 0; q < vect_index1; q++){ //start of for 19
            backdistrt[q] = 0;
          } //end of for 19
          //end of section corresponding to R code ******** if (x>1) {backdistrt[1:(x-1+N)] <- 0} 
          
          
        }       //end of top_level_if_4
        
        if (sum(backdistrt)==0) { //start of top_level_if_5
            dendistrt = dendistrt_N;      //correspoding to R code: dendistrt <- numeric(N) 
          }else{ //end of top_level_if_5, start of top_level_else_5
          
          backdistrt = backdistrt/sum(backdistrt);
          
          
          //start of section concerning R code ******** dendistrto <- dendistrto + c(rep(1,N),nichedistrt) ...
          
          w4 = 0;
          
          //also include nichdistrt in the same vector
          for(int y = N; y < vect_size; y++){ //start of for 6
            rep_vec2[y] = nichedistrt[w4];
            w4 = w4 + 1;
          } //end of for 6
          
          dendistrto = dendistrto + rep_vec2 * backdistrt * d * events_time[0];  
          
          //end of section concerning R code ******** dendistrto <- dendistrto + c(rep(1,N),nichedistrt) ...
          
          
          if (sum(dendistrto)==0) {  //start of top_level_if_6
               dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
            } //end of top_level_if_6
          
          dendistrto <- dendistrto/sum(dendistrto);
          
          
          
          
          
          
          
          //************************** start of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
          temp1 = N + 1;
          temp2 = 2 * N;
          int x_1 = 0;
          for(int p = temp1; p <= temp2; p++){ //start of for 3
            dendistrt[x_1] = dendistrto[p];
            x_1 = x_1 + 1;
          }  // end of for 3
          
          //************************** end of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
          
          
          //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
          Rcpp::NumericVector dendTemp2(N);
          for(int l = 0; l <= N; l++){ //start of for 4
            dendTemp2[l] = dendistrto[l];
          }  // end of for 4
          dendistrt[0] = dendistrt[0] + sum(dendTemp2);
          //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
          
          
          
          
          
          
        }       //end of top_level_else_5
        
      }  //end of top_level_else_3 
      
      
      //****************************************************************************************************************  
      
      
      
      
      
      
      
    }// end of top_level_if_2
    
    int tmp = -1; // deliberating declared within this if. This version of tmp exists only within outer if_1
    
    
    for(int i = 0; i < nrow_ev; i++){ //start of for loop
      
      dendistrt2 = dendistrt;
      
      
      
      if(events_type[i] == 1){ // start of outer_if_1
        
        events_n = events["n"];
        pt = events_n[i];
        
       
        //start of section to calculate cumsum
        for(int j = 0; j < dendistrto.size(); j++){ //start of for 2
          
          accum += dendistrto[j];
          tmp1[j] = accum;
          
          //end of section to calculate cumsum
          
          //concerning R code:  tmp <- min(which(tmp1>=abs(pt))) 
          if(tmp1[j] >= abs(pt) && tmp < 0){ 
            tmp = j;
          } //end of section concerning R code:  tmp <- min(which(tmp1>=abs(pt)))
          
          
        } //end of for 2
        
        accum = 0;       // to reset the cumsum
        
        //end of section to calculate cumsum
        
        
        if(pt<0){ //start of if_2
          
          
          //******************************this section concerns the R code if (tmp<(2*N)) {dendistrto[(tmp+1):(2*N)] <- 0}
          if(tmp<(2*N)){ //start of if_3
            
            temp1 = tmp + 1;
            temp2 = 2 * N;
            for(int k = temp1; k <= temp2; k++){ //start of for 2
              dendistrto[k] = 0;
            }  // end of for 2
            
            
          }   // *******************end of if_3 ****** this section concerns the R code if (tmp<(2*N)) {dendistrto[(tmp+1):(2*N)] <- 0}
          
          dendistrto[tmp] = dendistrto[tmp]-tmp1[tmp]+abs(pt); 
          
          
        }else{  //end of if_2 and start of if_2_else,i.e. the else corresponding to if_2
          
           if(tmp > 1){ //start of if_4: section concerning the R code    if (tmp>1) {dendistrto[1:(tmp-1)] <- 0}
            temp2 = tmp - 1;
            
            for(int k = 1; k <= temp2; k++){ //start of for 3
              dendistrto[k] = 0;
            } //end of for 3
            
          } //end of if_4  ****** end of if_4: section concerning the R code    if (tmp>1) {dendistrto[1:(tmp-1)] <- 0}
          
          
          dendistrto[tmp] = tmp1[tmp] - abs(pt); 
          
          
        } // end of if_2_else
        
       
        if (sum(dendistrto)==0 || sum(dendistrto<0)>0) { //start of if_5 
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break; 
        } //end of if_5
        
        dendistrto = dendistrto/sum(dendistrto);
        
        
        
        //************************** start of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        temp1 = N + 1;
        temp2 = 2 * N;
        int x1 = 0;
        for(int p = temp1; p <= temp2; p++){ //start of for 3
          dendistrt[x1] = dendistrto[p];
          x1 = x1 + 1;
        }  // end of for 3
        
        //************************** end of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        
        
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        Rcpp::NumericVector dendTemp(N);
        for(int l = 0; l <= N; l++){ //start of for 4
          dendTemp[l] = dendistrto[l];
        }  // end of for 4
        dendistrt[0] = dendistrt[0] + sum(dendTemp);
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        
        if (sum(dendistrt==dendistrt2) == N) { // start of if_6
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_6
        
        
        //start of section concerning R code ******** tmpo <- dendistrto*c(rep(1,N),)
       
        
        //Include nichdistrt in the same vector which already has rep(1,N)
        for(int y = N; y < vect_size; y++){ //start of for 6
          
          rep_vec2[y] = nichedistrt[w4];
          w4 = w4 + 1;
        } //end of for 6
        
        tmpo = dendistrto * rep_vec2;  
        
        //end of section concerning R code ******** tmpo <- dendistrto*c(rep(1,N),nichedistrt)
        
        
        
        
        //************************** start of section corresponding to tmpt <- tmpo[(N+1):(2*N)]
        temp1 = N + 1;
        temp2 = 2 * N;
        Rcpp::IntegerVector tmpt;
        x = 0;
        for(int q = temp1; q <= temp2; q++){ //start of for 7
          tmpt[x] = tmpo[q];
          x = x + 1;
        }  // end of for 7
        
        //************************** end of section corresponding to tmpt <- tmpo[(N+1):(2*N)]
        
        
        
        //section corresponding to R code ******** tmpt[1] <- tmpt[1]+sum(tmpo[1:N])
        Rcpp::IntegerVector tmpTemp(N);
        for(int l = 0; l <= N; l++){ //start of for 8
          tmpTemp[l] = tmpo[l];
        }  // end of for 8
        tmpt[0] = tmpt[0] + sum(tmpTemp);
        //section corresponding to R code ******** tmpt[1] <- tmpt[1]+sum(tmpo[1:N])
        
        
        
        //start of section corresponding to R code **** a <- a+events$a[z]*(sum(dendistrt*c(1:N))-sum(tmpt*c(1:N)))
        
        Rcpp::IntegerVector temp_vec3;
        Rcpp::NumericVector temp_vec4;
        Rcpp::IntegerVector temp_vec5;
        
        for(int t = 0; t < N; t++){ //start of for 8
          temp_vec3[t] = t + 1;
          
        } //end of for 8
        
        
        Rcpp::NumericVector double_temp_vec3 = Rcpp::NumericVector (temp_vec3); //typecast from IntegerVec to NumericVec    **********
        temp_vec4 = dendistrt * double_temp_vec3;
        temp_vec5 = tmpt * temp_vec3;
        
        a = a + events_a[i] * (sum(temp_vec4) - sum(temp_vec5));   
        
        //end of section corresponding to R code **** a <- a+events$a[z]*(sum(dendistrt*c(1:N))-sum(tmpt*c(1:N))) 
        
        
        
        //start of section corresponding to R code **** b <- b+events$b[z]*(sum(tmpt*c(1:N))-sum(dendistrt*c(1:N))) 
        
        b = b + events_b[i] * (sum(temp_vec5) - sum(temp_vec4));    
        
        //end of section corresponding to R code **** b <- b+events$b[z]*(sum(tmpt*c(1:N))-sum(dendistrt*c(1:N))) 
        
        
        // start of section corresponding to R code  nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
        
        for(int s = 0; s < N - 1; s++){
          nichedistrt[s] = pow(-1 * (s/exp(a) + pow(exp(1),-6)),(exp(b) + pow(exp(1),-6))); //****************************************
          
        }
        
        // end of section corresponding to R code  nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
        
        
        
        x = events_xmin[i];
        
        //equivalent to R code **** min(which(dendistrt>0))
        int temp_index = -1;
        for(int e = 0; e < dendistrt.size(); e++){ //start of for 9
          if(dendistrt[e] > 0 && temp_index < 0){
            temp_index = e;  // minimum index for which dendistrt > 0 will be stored here
          } 
        } //end of for 9
        
        if(x > temp_index) { //start of if_7
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        }   //end of if_7
        
        if(x > 1){ // start of if_8
          int temp_index2 = 0;
          temp_index2 = x-1+N;
          for(int f = 0; f<temp_index2;f++){ //start of for 10
            backdistrt[f] = 0;
          }//end of for 10
          
        } //end of if_8
        
        if(sum(backdistrt)==0){   //start of if_9
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_9
        
        backdistrt = backdistrt/sum(backdistrt);
        
        
        
        //start of section concerning R code ******** tmpo <- dendistrto*c(rep(1,N),nichedistrt)
        w4 = 0;
        
        //also include nichdistrt in the same vector
        for(int y = N; y < vect_size; y++){ //start of for 6
          
          rep_vec2[y] = nichedistrt[w4];
          w4 = w4 + 1;
        } //end of for 6
        
        Rcpp::NumericVector events_type_double = Rcpp::NumericVector(events_type);//typecasting as type double
        
        
        //start of section corresponding to R code **** dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
        if(i + 1 <= nrow_ev){
          dendistrto = dendistrto + rep_vec2 * backdistrt *  d * events_time[i+1] - events_type_double[i];
        }else{
          dendistrto = dendistrto + rep_vec2 * backdistrt *  d * edget - events_type[i];
        }
        //end of section corresponding to R code **** dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
        
        
        //end of section concerning R code ******** tmpo <- dendistrto*c(rep(1,N),nichedistrt)
        
        
        
        if (sum(dendistrto)==0 || sum(dendistrto<0)>0) { //test that this runs the way it is expected to***** start of if_10
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_10
        
        dendistrto <- dendistrto/sum(dendistrto);
        
        
        
        //************************** start of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        temp1 = N + 1;
        temp2 = 2 * N;
        int y1 = 0;
        for(int p = temp1; p <= temp2; p++){ //start of for 11
          dendistrt[y1] = dendistrto[p];
          y1 = y1 + 1;
        }  // end of for 11
        
        //************************** end of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        
        
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        Rcpp::IntegerVector dendTemp2(N);
        for(int l = 0; l <= N; l++){ //start of for 12
          dendTemp2[l] = dendistrto[l];
        }  // end of for 12
        dendistrt[0] = dendistrt[0] + sum(dendTemp2);
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        
      } //end of outer_if_1
      
      
      
      //*******************************************************************************************************************************************************************  
      
      
      if(events_type[i] == 2){ // start of outer_if_2
        
        int na = events_n[i];
        int na2 = 0;
        
        //start of section corresponding to R code ****** tmp <- numeric(2*N) 
        int temp_N = 2 * N;    //also consider declaring this at the start
        Rcpp::IntegerVector tmp2;
        int tmp3 = 0;
        for(int s = 0; s < temp_N; s++){  //start of for 13
          tmp2[s] = 0;    //i used tmp2 instead of tmp again as was used in the original R code. Did so to avoid confusion.
        }    //end of for 13
        //end of section corresponding to R code ****** tmp <- numeric(2*N)  
        
        if(na>0){     //start of if_12
          
          
          //start of section corresponding to ********* tmp[(na+1):(2*N)] <- dendistrto[1:(2*N-na)]
          temp1 = na + 1;
          temp2 = 2 * N; 
          temp3 = 1;
          temp4 = 2 * N - na;
          
          for(int m = 0; m < temp2; m++){ //start of for 14
            
            if(m<=temp4){ //only assign up to 2 * N - na
              tmp2[temp1]= dendistrto[m]; 
              temp1++;
            }
            if(m>temp4){
              temp_vec5 = dendistrto[m]; 
            }
          } //end of for 14
          //end of section corresponding to ********* tmp[(na+1):(2*N)] <- dendistrto[1:(2*N-na)] 
          
          tmp2[temp2] = tmp2[temp2] + sum(temp_vec5);
          
          
          
          //start of section concerning R code ******** dendistrto <- tmp*c(rep(1,N),nichedistrt)
          w4 = 0;
          
          //also include nichdistrt in the same vector
          for(int y = N; y < vect_size; y++){ //start of for 6
            
            rep_vec2[y] = nichedistrt[w4];
            w4 = w4 + 1;
          } //end of for 6
          
          Rcpp::NumericVector tmp2_double = Rcpp::NumericVector (tmp2); //typecasting
          dendistrto = tmp2_double * rep_vec2;
          
          //end of section concerning R code ******** dendistrto <- tmp*c(rep(1,N),nichedistrt)
          
        } //end of if_12
        
        if(na < 0){ //start of if_13
          
          na2 = abs(na);
          
          //start of section corresponding to ********* tmp[1:(2*N+na)] <- dendistrto[(-na+1):(2*N)] 
          temp1 = na2 + 1;
          temp2 = 2 * N; 
          temp3 = 0;
          temp4 = 2 * N + na;
          
          for(int m = temp1; m < temp2; m++){ //start of for 14
            
            if(m<=temp4){ //only assign up to 2 * N + na
              tmp2[temp3]= dendistrto[m]; 
              temp3++;
            }
            if(m>temp4){
              temp_vec5 = dendistrto[m]; 
            }
          } //end of for 14
          //end of section corresponding to ********* tmp[1:(2*N+na)] <- dendistrto[(-na+1):(2*N)] 
          
          
          //start of section corresponding to ********* tmp[1] <- tmp[1]+sum(dendistrto[1:(-na)])
          for(int g = 0; g < na2; g++){ // start of for 15
            temp_vec5 = dendistrto[g];
          }// end of for 15
          
          tmp2[1] = tmp2[1] + sum(temp_vec5);
          
          //end of section corresponding to ********* tmp[1] <- tmp[1]+sum(dendistrto[1:(-na)])
          
          dendistrto = tmp2;
          
        } //end of if_13
        
        
        if (sum(dendistrto)==0 || sum(dendistrto<0)>0) { //test that this runs the way it is expected to***** start of if_14
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_14
        
        dendistrto = dendistrto/sum(dendistrto);
        
        
        
        //************************** start of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        temp1 = N + 1;
        temp2 = 2 * N;
        int y1 = 0;
        for(int p = temp1; p <= temp2; p++){ //start of for 16
          dendistrt[y1] = dendistrto[p];
          y1 = y1 + 1;
        }  // end of for 16
        
        //************************** end of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        
        
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        Rcpp::IntegerVector dendTemp2(N);
        for(int l = 0; l <= N; l++){ //start of for 17
          dendTemp2[l] = dendistrto[l];
        }  // end of for 17
        dendistrt[0] = dendistrt[0] + sum(dendTemp2);
        //end of section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        
        u = u + na;
        Rcpp::IntegerVector temp_vec6;
        
        
        //start of section corresponding to tmp <- sum(dendistrt2*c(1:N))+na-sum(dendistrt*c(1:N))
        for(int h = 0; h < N; h++){
          temp_vec6[h] = h + 1;
        }
        Rcpp::NumericVector temp_vec6_double = Rcpp::NumericVector (temp_vec6); //typecasting
        tmp3 = sum(dendistrt2 * temp_vec6_double) + na - sum(dendistrt * temp_vec6_double);
        //end of section corresponding to tmp <- sum(dendistrt2*c(1:N))+na-sum(dendistrt*c(1:N))
        
        if(na > 0 && tmp3 > 0){  //start of if_15
          a = a + events_a[i] * tmp3;
          b = b + events_b[i] * tmp3;
          
          
          
          // start of section corresponding to R code  nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
          
          for(int s = 0; s < N - 1; s++){
            nichedistrt[s] = pow(-1 * (s/exp(a) + pow(exp(1),-6)),(exp(b) + pow(exp(1),-6))); //****************************************
            
          }
          
          // end of section corresponding to R code  nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
          
          
        }  //end of if_15
        
        // start of section corresponding to backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
        
        backdistrt = pnorm(Rcpp::seq(-1 * (N - 1),N),u,sd) - pnorm(Rcpp::seq(-1 * N,N - 1),u,sd);
        
        // End of section corresponding to backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
        
        
        x = events_xmin[i];
        
        // start of section equivalent to R code **** min(which(dendistrt>0)) **********************
        int temp_index = -1;
        for(int e = 0; e < dendistrt.size(); e++){ //start of for 18
          if(dendistrt[e] > 0 && temp_index < 0){
            temp_index = e;  // minimum index for which dendistrt > 0 will be stored here
          } 
        } //end of for 18
        // end of section equivalent to R code **** min(which(dendistrt>0)) ************************
        
        if(x > temp_index){ //start of if_16
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_16
        
        if(x > 1){ //start of if_17
          
          //start of section corresponding to R code ******** if (x>1) {backdistrt[1:(x-1+N)] <- 0} 
          int vec_index1 = x - 1 + N;
          for(int q = 0; q < vec_index1; q++){ //start of for 19
            backdistrt[q] = 0;
          } //end of for 19
          //end of section corresponding to R code ******** if (x>1) {backdistrt[1:(x-1+N)] <- 0} 
          
        }//end of if_17
        
        if(sum(backdistrt)==0){ //start of if_18
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
          
        }  //end of if_18
        
        backdistrt = backdistrt/sum(backdistrt);
        
        
        //start of section concerning R code ******** tmpo <- dendistrto*c(rep(1,N),nichedistrt)
         w4 = 0;
        
        //also include nichdistrt in the same vector
        for(int y = N; y < vect_size; y++){ //start of for 6
          
          rep_vec2[y] = nichedistrt[w4];
          w4 = w4 + 1;
        } //end of for 6
        
  
        //start of section corresponding to R code **** dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
        if(i + 1 <= nrow_ev){
          dendistrto = dendistrto + rep_vec2 * backdistrt *  d * events_time[i+1] - events_type[i];
        }else{
          dendistrto = dendistrto + rep_vec2 * backdistrt *  d * edget - events_type[i];
        }
        //end of section corresponding to R code **** dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
        
        
        if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {  //start of if_19
          dendistrt = dendistrt_N;      //corresponding to R code: dendistrt <- numeric(N) 
          break;
        } //end of if_19
        
        
        dendistrto = dendistrto/sum(dendistrto);
        
        //************************** start of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        temp1 = N + 1;
        temp2 = 2 * N;
        int y2 = 0;
        for(int p = temp1; p <= temp2; p++){ //start of for 20
          dendistrt[y2] = dendistrto[p];
          y2 = y2 + 1;
        }  // end of for 20
        
        //************************** end of section corresponding to dendistrt <- dendistrto[(N+1):(2*N)]
        
        
        //section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        Rcpp::IntegerVector dendTemp3(N);
        for(int l = 0; l <= N; l++){ //start of for 21
          dendTemp3[l] = dendistrto[l];
        }  // end of for 21
        dendistrt[0] = dendistrt[0] + sum(dendTemp3);
        //end of section corresponding to R code ******** dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
        
        
      }//end of outer_if_2
      
      
      
      
    } //end of for loop
    
    
  }  // end of top_level_if_1
  
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("dendistrt") = dendistrt , Rcpp::_["dendistrto"] = dendistrto, 
                                         Rcpp::_["nichedistrt"] = nichedistrt,Rcpp::_["backdistrt"] = backdistrt,Rcpp::_["a"] = a,Rcpp::_["b"] = b,
                                                 Rcpp::_["u"] = u,Rcpp::_["x"] = x);
  
 return result;
  
  
}

