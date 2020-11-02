# Gennect2020
Local ver. 으로 수정됨. 

실행방법 

0. git clone으로 전체 코드 다운로드 
1. pip install -r requirements.txt 로 필요한 패키지 설치 
2. app 폴더로 이동 후 python main.py 로 서버 실행 
3. localhost:5000 에 접속하여 확인 

%%% 주의사항 

현재 개발 버젼으로, 대소문자 구별하여 hgnc_symbol을 넣지 않으면 에러가 뜸. 

예를 들어, "BATF" 는 가능하나 "batf"는 불가능하며, "TP53"은 가능하나 "p53"은 불가능함. 
