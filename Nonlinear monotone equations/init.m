% init 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [nprob,n,m,x0]=initf(nprob)
% This function sets n,m, and the standard starting    
% point based on the nprob and returns it to initpt     
% function.                                                                                                    
% Created on 10/30/94 by Madhu Lamba                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nprob,n,x0] = init(NO)

%global FIRSTIME

switch NO
    %% nprob='problem1'; cuter
      case 1
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 2
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 3
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 4
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem1'; cuter
      case 5
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem1'; cuter
      case 6
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem1'; cuter
      case 7
        nprob='problem1';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
    
    %% nprob='problem1'; cuter
      case 8
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 9
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 10
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 11
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem1'; cuter
      case 12
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem1'; cuter
      case 13
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem1'; cuter
      case 14
        nprob='problem1';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
        
    %% nprob='problem1'; cuter
      case 15
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 16
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 17
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem1'; cuter
      case 18
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem1'; cuter
      case 19
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem1'; cuter
      case 20
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem1'; cuter
      case 21
        nprob='problem1';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
        
    %% nprob='problem2'; cuter
      case 22
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem2'; cuter
      case 23
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 24
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 25
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem2'; cuter
      case 26
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem2'; cuter
      case 27
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem2'; cuter
      case 28
        nprob='problem2';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem2'; cuter
      case 29
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem2'; cuter
      case 30
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 31
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 32
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem2'; cuter
      case 33
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem2'; cuter
      case 34
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem2'; cuter
      case 35
        nprob='problem2';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem2'; cuter
      case 36
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem2'; cuter
      case 37
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 38
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem2'; cuter
      case 39
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem2'; cuter
      case 40
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem2'; cuter
      case 41
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem2'; cuter
      case 42
        nprob='problem2';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 
        
    %% nprob='problem3'; cuter
      case 43
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem3'; cuter
      case 44
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 45
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 46
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem3'; cuter
      case 47
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem3'; cuter
      case 48
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem3'; cuter
      case 49
        nprob='problem3';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end         
        
%% nprob='problem3'; cuter
      case 50
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem3'; cuter
      case 51
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 52
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 53
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem3'; cuter
      case 54
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem3'; cuter
      case 55
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem3'; cuter
      case 56
        nprob='problem3';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end        
 
    %% nprob='problem3'; cuter
      case 57
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem3'; cuter
      case 58
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 59
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem3'; cuter
      case 60
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem3'; cuter
      case 61
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem3'; cuter
      case 62
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem3'; cuter
      case 63
        nprob='problem3';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end

    %% nprob='problem4'; cuter
      case 64
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem4'; cuter
      case 65
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 66
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 67
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem4'; cuter
      case 68
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem4'; cuter
      case 69
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem4'; cuter
      case 70
        nprob='problem4';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem4'; cuter
      case 71
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem4'; cuter
      case 72
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 73
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 74
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem4'; cuter
      case 75
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem4'; cuter
      case 76
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem4'; cuter
      case 77
        nprob='problem4';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem4'; cuter
      case 78
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem4'; cuter
      case 79
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 80
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem4'; cuter
      case 81
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem4'; cuter
      case 82
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem4'; cuter
      case 83
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem4'; cuter
      case 84
        nprob='problem4';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end         
        
     %% nprob='problem5'; cuter
      case 85
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem5'; cuter
      case 86
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 87
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 88
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem5'; cuter
      case 89
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem5'; cuter
      case 90
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem5'; cuter
      case 91
        nprob='problem5';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem5'; cuter
      case 92
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem5'; cuter
      case 93
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 94
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 95
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem5'; cuter
      case 96
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem5'; cuter
      case 97
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem5'; cuter
      case 98
        nprob='problem5';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem5'; cuter
      case 99
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem5'; cuter
      case 100
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 101
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem5'; cuter
      case 102
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem5'; cuter
      case 103
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem5'; cuter
      case 104
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem5'; cuter
      case 105
        nprob='problem5';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end        
        
     %% nprob='problem6'; cuter
      case 106
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem6'; cuter
      case 107
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 108
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 109
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem6'; cuter
      case 110
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem6'; cuter
      case 111
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem6'; cuter
      case 112
        nprob='problem6';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem6'; cuter
      case 113
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem6'; cuter
      case 114
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 115
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 116
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem6'; cuter
      case 117
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem6'; cuter
      case 118
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem6'; cuter
      case 119
        nprob='problem6';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem6'; cuter
      case 120
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem6'; cuter
      case 121
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 122
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem6'; cuter
      case 123
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem6'; cuter
      case 124
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem6'; cuter
      case 125
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem6'; cuter
      case 126
        nprob='problem6';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end  

%% nprob='problem7'; cuter
      case 127%64%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem7'; cuter
      case 128%65%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 129%66%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 130%67%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem7'; cuter
      case 131%68%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem7'; cuter
      case 132%69%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem7'; cuter
      case 133%70%
        nprob='problem7';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem7'; cuter
      case 134%71%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem7'; cuter
      case 135%72%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 136%73%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 137%74%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem7'; cuter
      case 138%75%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem7'; cuter
      case 139%76%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem7'; cuter
      case 140%77%
        nprob='problem7';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem7'; cuter
      case 141%78%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem7'; cuter
      case 142%79%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 143%80%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem7'; cuter
      case 144%81%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem7'; cuter
      case 145%82%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem7'; cuter
      case 146%83%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem7'; cuter
      case 147%84%
        nprob='problem7';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
        
    %% nprob='problem8'; cuter
      case 148
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem8'; cuter
      case 149
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 150
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 151
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem8'; cuter
      case 152
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem8'; cuter
      case 153
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem8'; cuter
      case 154
        nprob='problem8';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem8'; cuter
      case 155
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem8'; cuter
      case 156
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 157
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 158
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem8'; cuter
      case 159
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem8'; cuter
      case 160
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem8'; cuter
      case 161
        nprob='problem8';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem8'; cuter
      case 162
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem8'; cuter
      case 163
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        x0(1)=1/3;
        for i=1:n-1
            x0(i+1)=1/3*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 164
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem8'; cuter
      case 165
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem8'; cuter
      case 166
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem8'; cuter
      case 167
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem8'; cuter
      case 168
        nprob='problem8';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
        
%     %% nprob='problem10'; cuter
%       case 127
%         nprob='problem10';
%         n=4;
%         x0=zeros(n,1);   
%         
%     %% nprob='problem10'; cuter
%       case 128
%         nprob='problem10';
%         n=4;
%         x0=[3;0;0;0];  
%         
%     %% nprob='problem10'; cuter
%       case 129
%         nprob='problem10';
%         n=4;
%         x0=[1;1;1;0]; 
%         
%     %% nprob='problem10'; cuter
%       case 130
%         nprob='problem10';
%         n=4;
%         x0=[0;1;1;1]; 

    %% nprob='problem11'; cuter
      case 169
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem11'; cuter
      case 170
        nprob='problem9';
        n=10000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem11'; cuter
      case 171
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem11'; cuter
      case 172
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem11'; cuter
      case 173
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem11'; cuter
      case 174
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem11'; cuter
      case 175
        nprob='problem9';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem11'; cuter
      case 176
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem11'; cuter
      case 177
        nprob='problem9';
        n=50000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem11'; cuter
      case 178
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem11'; cuter
      case 179
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem11'; cuter
      case 180
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem11'; cuter
      case 181
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem11'; cuter
      case 182
        nprob='problem9';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem11'; cuter
      case 183
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem11'; cuter
      case 184
        nprob='problem9';
        n=100000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem11'; cuter
      case 185
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem11'; cuter
      case 186
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem11'; cuter
      case 187
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem11'; cuter
      case 188
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem11'; cuter
      case 189
        nprob='problem9';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
        
%% nprob='problem12'; cuter
      case 190%85%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
    
    %% nprob='problem12'; cuter
      case 191%86%
        nprob='problem12';
        n=10000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem12'; cuter
      case 192%87%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem12'; cuter
      case 193%88%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem12'; cuter
      case 194%89%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem12'; cuter
      case 195%90%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem12'; cuter
      case 196%91%
        nprob='problem12';
        n=10000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end    
    
    %% nprob='problem12'; cuter
      case 197%92%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
    
    %% nprob='problem12'; cuter
      case 198%93%
        nprob='problem12';
        n=50000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem12'; cuter
      case 199%94%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem12'; cuter
      case 200%95%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem12'; cuter
      case 201%96%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem12'; cuter
      case 202%97%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem12'; cuter
      case 203%98%
        nprob='problem12';
        n=50000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end 

    %% nprob='problem12'; cuter
      case 204%99%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
    
    %% nprob='problem12'; cuter
      case 205%100%
        nprob='problem12';
        n=100000;
        x0=0.1*ones(n,1);
        
    %% nprob='problem12'; cuter
      case 206%101%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
        x0(1)=0.5;
        for i=1:n-1
            x0(i+1)=0.5*x0(i);
        end
        
    %% nprob='problem12'; cuter
      case 207%102%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=(i-1)/n;
        end
    
    %% nprob='problem12'; cuter
      case 208%103%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1/i;
        end
        
    %% nprob='problem12'; cuter
      case 209%104%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=i/n;
        end
        
    %% nprob='problem12'; cuter
      case 210%105%
        nprob='problem12';
        n=100000;
        x0=ones(n,1);
        for i=1:n
            x0(i)=1-i/n;
        end
end