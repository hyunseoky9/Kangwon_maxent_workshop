# 이번에는 민물에 사는 종을 maxent 를 통해 그려보도록 하자.
# 육상종과 달리 민물에 사는 종들은 하천에 주로 서식하기 때문에 면이 아닌 선형적인
# 서식공간에 분포해 있다. 그렇기 때문에 주로 면적을 다루는 레스터 형식의 예측변수
# 데이터는 적합하지 않다. 그렇다 하더라도 데이터 형식만 다른 것이지 Maxent 의 
# 원리는 종의 서식환경에 관계없이 똑같이 적용된다.

# 이 예시를 통해서 하천 네트워크에 서식하는 종들의 데이터를 어떻게 Maxent 를 실행하는
# 함수에 넣고 분포도를 시각화하는지 알아보자. 대부분의 절차는 전 예시와 같으므로
# 간단하게 짚고 넘어가고 하천 네트워크를 이용해야해서 다른점만 자세히 다루고자 한다.

# 또 이번에는 미래 기후에 따라 변한 환경을 반영한 예측변수들을 이용해
# 현재 환경을 반영한 데이터로 만들어진 Maxent 모델로 미래의 분포도를 예측하는 연습도 해보자.
# 미래 기후변화에 따른 종들의 미래 서식지를 예측하는 것은 ENM의 대표적인 활용방법중
# 하나이다.

# 이번에 예시로 쓰일 종은 제첩 (Corbicula fluminea)으로 민물에 사는 패류이며 
# 미국에서는 외래종이다.
# 현제 미국의 Tennessee water resource region (두자리 Hydrological Unit Code 에서 6에 해당하는 지역)
# 에서 제첩의 분포도와 미래의 분포도를 예측해보자.


# 목차: 
# 1. 데이터 불러오기
# 2. ENMeval을 통한 Maxent 실행 (Running Maxent using ENMEval)
# 3. 모델 선택
# 4. 분포도 예측

# 1. 발생 데이터
# 먼저 현재 지정되있는 매개들을 전부다 없애자. 
rm(list=ls())
# Rstudio를 사용하고 있으면 Environment 탭이 지워진 것을 볼 수 있다.
# 항상 새로운 작업을 시작하면 전에 사용하던 매개들을 없애고 시작하는 것이 좋은
# 습관이다. 나중에 혹시라도 전에 사용하던 매개들이 쓰이게 되면 에러가 날 경우가 
# 있기때문이다.

# 패키지들을 불러온다. 아직 지난 예시에서 받지 않은 패키지들이 있는데
# 다운받은후 불러오자.
#devtools::install_github("johnbaums/rmaxent")
#install.packages('data.table')
#install.packages('usdm')
library(ENMeval)
library(sf)
library(rmaxent) # 이번엔 ENMevaluate 을 쓸때 maxnet 대신 maxent.jar 을 써보기로 하자.
library(data.table)
library(dplyr)
library(tidyverse)
library(usdm)
library(rstudioapi)
# 작업 디랙토리를 셋팅한다.
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))


# 종의 발생데이터를 불러온다.
spdata = fread('./data/Corbicula fluminea_occ_current.csv')

# 데이터를 잠깐 보자
head(spdata)
# 아까 예시보다 훨씬더 행들이 많은 것을 볼 수 있다. 종의 발생 지역좌표 (decimalLongitude/Latitude)
# 만 있는 것이 아니라 그 좌표에 해당하는 하천 네트워크 아이디 (comid),
# 기록일, GBIF ID, 좌표가 있는 HUC (huc12,huc2) 코드 등의 데이터 한열마다 기별할수 있는 정보들이 있다.
# 실제 연구에서 쓰이는 데이터들은 이와같이 좌표만 있지않고 여러 정보들이 더 있는 경우가 많다.
# 그 외의 행들은 발생 지역에 해당하는 환경 예측변수들의 값이다. 
# 전 예시에서는 발생 데이터와 배경 데이터는 좌표만 있고 환경 데이터는 따로 있었지만
# 데이터가 방대한 경우에는 환경데이터와 발생데이터를 매칭시키는 과정이 오래 걸릴 수 있
# 기 때문에 좌표데이터와 같이 미리 환경데이터도 붙여놓을 수도 있다.

# 다음은 연구지역에 있있는 데이터를 불러오자.
bgdata = fread('./data/Corbicula fluminea_projarea_current.csv')
head(bgdata)

# 연구 지역은 모든 발생데이터를 포괄할 수 있는 최소한으로 구성된 6자리 HUC 유닛들로한다.
# 위 데이터는 그 HUC 유닛들 안에 있는 모든 하천 유닛들을 포함한다. 하천유닛들은 
# National Hydrography Dataset에의해 정의되어 있는 유닛들이며 각 유닛마다 ID가 
# 존재한다 (comid).


# 현재 지하수 유입율 (BFI) 과 하천 순서 (stream order) 예측변수 들중에는
# 각각 누락된 값이 있다. 이들을 데이터에서 제외하자.
if(length(which(is.na(spdata$BFI)))>0)
{
  spdata = spdata[-which(is.na(spdata$BFI)),]
}
if(length(which(is.na(bgdata$BFI)))>0)
{
  bgdata = bgdata[-which(is.na(bgdata$BFI)),]    
}
# exclude rows with streamorder value of -9
if(length(which(spdata$streamorder==-9))>0)
{
  spdata = spdata[-which(spdata$streamorder==-9),]  
}
if(length(which(bgdata$streamorder==-9))>0)
{
  bgdata = bgdata[-which(bgdata$streamorder==-9),]
}

# 데이터를 Tennessee water resource region (HUC 06)에 그래프해보자.
# 우선 지역 공간정보를 담고 있는 파일을 부른다.
huc06 = st_read('./data/huc_06.shp')
huc06 = st_transform(huc06,4326)

# 발생 데이터에서 좌표를 huc06와 같은 crs를 가진 공간 매개로 지정후 그린다.
spoccdatapts = st_as_sf(spdata, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
plot(spoccdatapts$geometry,xlim=c(min(spdata$decimalLongitude)-0.1,max(spdata$decimalLongitude)+0.1),col='white')
plot(huc06,col='white',main='',add=TRUE)
plot(spoccdatapts$geometry,add=TRUE)

# 연구 지역에 있는 모든 데이터도 그려보자
bgdatapts = st_as_sf(bgdata, coords = c('decimalLongitude','decimalLatitude'), crs = 4326)
plot(bgdatapts$geometry,add=TRUE)
# 연구지역이 꽉찬 데이터로 꽉찬 것을 볼 수 있다.

# 앞서 보앗듯 환경에 관련된 예측변수들이 여럿있다. 이중에는 상관관계가 높은 변수들도
# 있기 때문에 다중공선성 분석(multicolinearity analysis) 을 통해 그들을 걸러낸다.
# 상관관계가 높은 변수들이 함께 Maxent 모델에 들어가면 모델 계수측정이
# 불안정하게 이뤄지기 때문에 좋지 않다.

# VIF(분산 팽창 요인, Variance Inflation Factor) 분석은 다중 회귀 분석에서 독립
# 변수 간의 다중공선성을 평가하는 데 사용되는 기법이다. 
# VIF 값은 각 변수를 다른 모든 변수로 회귀분석했을 때 얻어지는 결정 계수(R-squared)를
# 기반으로 계산된다. VIF 값이 1에 가까울수록 해당 독립 변수는 다른 독립 변수와의
# 상관성이 낮다는 것을 의미하며, 10 이상의 높은 VIF 값은 해당 변수가 다른 변수와
# 높은 상관성을 가지고 있음을 나타냅니다.
# 높은 VIF 값을 가진 변수는 모델에서 제거하거나 다른 변수와의 조합을 다시 고려하여
# 다중공선성 문제를 해결할 수 있습니다. VIF 분석은 모델의 정확성과 해석 가능성을
# 높이기 위해 중요하게 여겨집니다.

usdm::vif(dplyr::select(bgdata,waterbody,BFI,streamorder,numday_above_tmax,numday_below_tmin,
           dd90_5c,dd90_8c,dd90_10c,dd120_5c,dd120_8c,dd120_10c,
           dd150_5c,dd150_8c,dd150_10c,avgtemp,maxflow,maxflowdate,
           minflow,minflowdate,avgflow))
# 보면 VIF 값이 10을 가볍게 넘는 변수들이 많다. 그만큼 상관관계가 큰 변수들이 많이 들어있다는 거다.

# 여러가지 조합들의 VIF 를 계산해본 후 저자는 다음과 같은 변수들의 조합으로 추려냈다.
usdm::vif(dplyr::select(bgdata,BFI,waterbody,streamorder,numday_above_tmax,dd90_8c,minflow,minflowdate,maxflowdate))

# 각각의 예측변수들의 설명은 다음과 같다.
# BFI: 지하수 유입률
# waterbody: 호수등 정적인 물에 있는 여부 (0,1 값)
# streamorder: 하천 순서
# numday_above_tmax: 연간 임계 온도 최대치 를 넘는 날들의 수 연평균.
# numday_below_tmin: 연간 임계 온도 최저치 를 넘는 날들의 수 연평균.
# dd90_8c: 쥴리안 날로 90일에 베이스 온도를 섭씨8도로 지정했을 때의 도일수 연평균
# minflow: 최저 유속 연평균
# minflowdate: 최저유속일 연평균
# maxflowdate: 최대유속일 연평균
# 위에서 온도와 관련된 예측 변수들은 시뮬레이션을 통하여 만든 하천 표면온도로 만들었다.
# 하천 생물은 하천의 온도에 영향을 받기때문에 대기온도보다 하천온도를 쓰는것이 정확하다.
# 또한, 그 종의 관점에서 환경에대한 해석을 담은 변수들을 사용하는 것이 좋다.
# 얘룰둘오 numday_above_tmax는 수온을 종의 임계 온도 최대치에 반영하여 해석한 변수이다.
# 이런 종의 관점에서 환경을 해석한 변수들이 흔히쓰는 연평균 수온보다 ENM을 하는데 
# 있어서 더 정확성을 부여할 수 있다.


# 이제 다중공선성 분석후 추려낸 예측변수들만 데이터에 남겨놓자 
bgdata = data.frame(dplyr::select(bgdata,decimalLongitude,decimalLatitude,BFI,waterbody,streamorder,numday_above_tmax,dd90_8c,minflow,minflowdate,maxflowdate))
occ = data.frame(dplyr::select(spdata,decimalLongitude,decimalLatitude,BFI,waterbody,streamorder,numday_above_tmax,dd90_8c,minflow,minflowdate,maxflowdate))

# 이제 연구지역에서 10,000 개의 포인트를 무작위로 뽑에 배경 데이터로 삼는다.
if(nrow(bgdata)<10000)
{
  samplesize = nrow(bgdata)
} else {
  samplesize = 10000  
}
bg = bgdata[sample(1:nrow(bgdata), samplesize),]

#waterbody 예측변수는 범주형 변수이므로 factor 로 선언한다.
bg$waterbody = as.factor(bg$waterbody)
occ$waterbody = as.factor(occ$waterbody)

# ENMEval이 좌표가 있는 행들을 알아보기 위해
# 좌표가 있는 행들의 이름을 'longitude' 과 'latitude'로 바꾼다.
names(bg)[which(names(bg)=='decimalLongitude')] = 'longitude'
names(occ)[which(names(occ)=='decimalLongitude')] = 'longitude'
names(bg)[which(names(bg)=='decimalLatitude')] = 'latitude'
names(occ)[which(names(occ)=='decimalLatitude')] = 'latitude'


# 2. ENMeval을 통한 Maxent 실행 (Running Maxent using ENMEval)
# 이제 ENMevaluate 함수를 통해 준비된 데이터로 Maxent를 실행하고 튜닝을 하자.
# 튜닝을 하기위해 여러 RM과 feature class 매개변수값들을 지정해준다.
tune.args <- list(fc = c("L", "LQ","LQH","LQHP","LQHPT"), rm = 1:5) # model tuning parameter range
# 밑에 줄은 ENMevaluate 으로 maxent를 실행시킨다.
# Maxent를 실행하기 위한 방법으로는 이번에는 maxent.jar을 써보자.
# 데이터 분할 방법으로는 블록방식을 쓰자.
# 또 병렬프로그래밍을 통해 더빨리 모델을 만들어낼수도 있지만 지금은 꺼두자 (parallel=FALSE)
# 환경 예측변수를 넣었던 파라미터가 없는것을 볼 수 있다.
# 하천 생물은 대게 서식처가 면으로 묘사를 하지 않고 선으로 묘사하기 때문에 레스터데이타로
# 예측변수 파일을 만들지 않는다. 
# 보통은 그냥 예측변수를 테이블형식으로 가지고 있는다 (tabular). 테이블 형식의
# 데이터로 maxent를 ENMevaluate 을 실행하려면 다음과 같은 양식을 쓴다.
#e = ENMevaluate(occ, bg = bg, algorithm = "maxent.jar", tune.args = tune.args, 
#                partitions = "block",parallel=FALSE)
# 위 함수를 실행 시키려면 컴퓨터 사양에 따라 시간이 오래 걸릴 수 있으므로 미리 돌려서
# 저장해 놓은 매개를 불러오자
e = read_rds('./data/Corbicula_fluminea.rds')

# 여러 통계값을 그래프해보자.
evalplot.stats(e = e, stats = c("or.mtp", "auc.val","cbi.val"), color = "fc", x.var = "rm",
               error.bars = FALSE)

# 3. 모델 선택

# 전 예시와 달리 이번에는 Continuous Boyce Index (CBI) 가 가장 높은 모델을 선택한다.
res <- eval.results(e)
res2 = res[which(!is.na(res$cbi.val.avg)),]
opt.seq <- res2 %>% 
  filter(cbi.val.avg == max(cbi.val.avg))

# 선택된 모델의 매개변수:
opt.seq[1:2]
# 선택된 모델을 다른 매게에 따로 저장한다.
mod.seq <- eval.models(e)[[opt.seq$tune.args]]

# maxent.jar 로 Maxent를 실행하면 여러 결과를 보여주는 html이 만들어진다.
# AUC-ROC 그래프와 예측변수마다 분포도를 설명하는데 있어서의 중요도를 볼 수 있다.
# 하천 순서와 도일수가 가장 중요한 변수인 것을 볼 수 있다.
mod.seq

# 4. 분포도 예측
# 이제 선택한 모델로 Tennessee water resource region 에서 종의 현제 잠재 분포도도를 
# 예측해 보자.

# 우선 이 지역에의 예측변수를 가져온다. 
proj = fread('./data/huc6_projarea_tempNflow_predictors_current.csv')

# 아까와 같이 누락된 값이 있는 열들은 없애준다.
if(length(which(is.na(proj$BFI)))>0)
{
  proj = proj[-which(is.na(proj$BFI)),]
}
if(length(which(proj$streamorder==-9))>0)
{
  proj = proj[-which(proj$streamorder==-9),]
}

# 모델을 만들때 썼던 예측변수만 가져간다.
proj = data.frame(dplyr::select(proj,comid,BFI,waterbody,streamorder,numday_above_tmax_30.000,dd90_8c,minflow,minflowdate,maxflowdate,avgflow))
# 모델을 만들때 썼던 예측변수의 네이밍과 똑같은 이름을 사용한다.
names(proj)[5] = 'numday_above_tmax'

# project 함수를 통해 모델 예측을 실행한다.
projpred = rmaxent::project(mod.seq,proj[,2:ncol(proj)])
# 3가지 예측값이 나온 것을 볼 수 있다. (raw, logistic, cloglog')
str(projpred)
# 이중에 하천단위마다 종이 서식할 확률을 나타내줄 수 있는 cloglog를 사용하자.
projpred = projpred$prediction_cloglog

# 이제 이 확률에 어떤 확률 위로는 이곳에 서식한다고 정의해야한다.
# 그 기준점을 정하는 데는 여러 방법이 있는데 이중 우리는 MaxSSS Threshold 기준점을
# 사용하자. MaxSSS threshold 는 training 데이터의
# True Positive Rate 과 True Negative Rate의 합을 극대화 시켜주는 기준점이다.
# 다음 줄들을 통해 MaxSSS 를 구하자
occpred = rmaxent::project(mod.seq,occ[,3:ncol(occ)])
occpred = occpred$prediction_cloglog
bgpred = rmaxent::project(mod.seq,bg[,3:ncol(bg)])
bgpred = bgpred$prediction_cloglog
Pred = c(occpred,bgpred)
Sp.occ = c(rep(1,length(occpred)),rep(0,length(bgpred)))
maxsss = ecospat::ecospat.max.tss(Pred, Sp.occ)
maxsss = maxsss$max.threshold #maxsss threshold
maxsss # MaxSSS 기준점 값.

# 이제 이 기준점보다 높은 서식 확률값을 가진 곳은 서식한다(1),
# 낮은 값을 가진 곳은 서식하지 않는다(0)로 정의하자.
s_maxsss = rep(0, length(projpred))  
s_maxsss[which(projpred>maxsss)] = 1

# 이제 이 분포도를 시각화 해보자.
# 우선 지역 하천네트워크의 공간 데이터를 불러온 다음 그려본다.
flowlines = st_read('./data/huc_06_nhdflowlines.shp')
plot(flowlines$geometry,col='black')


# 하천 유닛중 서식하는 곳으로 정의된 곳을 조금더 두껍고 다른색의 라인으로 나타낸다
plot(flowlines$geometry[which(s_maxsss==1)],col='red',lwd=2, add=TRUE)
# 동쪽으로 더 많이 분포해 있는것을 볼 수 있다.

# 이번에는 미래의 기후에서의 잠재적 분포도를 알아보자
# 다음은 시뮬레이션을 통해 얻은 미래의 기후데이터를 기반으로 만든 이 지역에서의 예측변수다. 
projf = fread('./data/huc6_projarea_tempNflow_predictors_future.csv')

# 똑같은 방법으로 데이터를 가공한다.
if(length(which(is.na(projf$BFI)))>0)
{
  projf = projf[-which(is.na(projf$BFI)),]
}
if(length(which(projf$streamorder==-9))>0)
{
  projf = projf[-which(projf$streamorder==-9),]
}
projf = data.frame(dplyr::select(projf,comid,BFI,waterbody,streamorder,numday_above_tmax_30.000,dd90_8c,minflow,minflowdate,maxflowdate,avgflow))
names(projf)[5] = 'numday_above_tmax'
projfpred = rmaxent::project(mod.seq,projf[,2:ncol(projf)])
projfpred = projfpred$prediction_cloglog

# 아까 구한 MaxSSS threshold 값으로 하천 유닛마다 서식의 유무를 정의한다.
s_maxsssf = rep(0, length(projfpred))  
s_maxsssf[which(projfpred>maxsss)] = 1
plot(flowlines$geometry[which(s_maxsssf==1)],col='blue',lwd=2, add=TRUE)
# 미래 잠재적 분포도는 훨씬더 넓은 것을 볼 수 있다.

# 마지막으로 현재 분포도를 빨간색으로, 미래에 추가된 분포지역을 초록색으로 나타내 보자.
# 먼저 연구지역 테두리를 부른다.
s_maxsss_fsubc = s_maxsssf - s_maxsss
s_maxsss_fsubc[which(s_maxsss_fsubc==-1)] = 0
plot(huc06$geometry)
plot(flowlines$geometry[which(s_maxsss_fsubc==1)],col='green',lwd=1, add=TRUE)
plot(flowlines$geometry[which(s_maxsss==1)],col='red',lwd=1, add=TRUE)


# 현재 미래 분포도 들을 .shp 파일로 저장해서 후에 QGIS나 ArcGIS 로 작업을 할 수도 있다.
st_write(flowlines[which(s_maxsss==1),],'./predicted_current_distribution.shp') # 현재 분포도
st_write(flowlines[which(s_maxsss_fsubc==1),],'./predicted_future_distribution.shp') # 미래 분포도도


