# *이 실습의 코드와 커멘트들은 ENMEval2.0 예시설명서에 있는 내용을 가져와 번역하고 조금 더 덧붙여
# 만든 것입니다. (소스: https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html)

# ENMEval2.0 패키지를 통해 Maxent 알고리즘으로 Ecological Niche Modeling (ENM) 을 해보자.
# 본래 MaxEnt 알고리즘 프로그램은 육상 동식물의 분포도를 예측하는 기반으로 만들어 졌기 때문에
# 서식환경에 대한 예측변수 데이터들이 레스터 형식인 경우가 많다.
# 이 스크립트에서는 세발가락 나무늘보 발생 데이터를 가지고 이 종의 분포도를 예측해 보자.
# ENM 을 하기위해서 중요한 과정들을 순서대로 하기위해 과정들을 챕터별로 나누었다. 

# 목차: 
# 1. 데이터 취득과 가공 (Data Acquisition & Pre-processing)
# 2. 모델 평가를 위한 발생 데이터 분할 (Partitioning Occurrences for Evaluation).
# 3. ENMeval을 통한 Maxent 실행 (Running Maxent using ENMEval)
# 4. 모델 선택 (Model selection)

# 필드에서 가져온 순수한 발생 데이터를 데이터베이스에서 가져오고 그 데이터를 Maxent알고리즘 패키지에
# 대입하기 위해 가공하는 과정을 pre-processing 이라고 한다.

# 우선 필요한 R 패키지들을 불러온다. 불러오는 순서가 중요하다. 한패키지의 함수가 다른 패키지의 동명의 함수들을
# 덮어씌울수도 있기 때문이다.
# 패키지를 아직 다운로드 받지 않았다면 아래 install.package가 있는 줄들을 uncomment하고 먼저 받는다.
# install.packages('ENMeval')
# install.packages('raster')
# install.packages('tidyverse')
# install.packages('dplyr')

library(ENMeval)
library(raster)
library(dplyr)
library(tidyverse)

# 다시 똑같은 분석을 할 수 있로록 랜덤함수 고정값을 주자.
set.seed(48)

# 세발가락 나무늘보 (Bradypus variegatus)의 발생 (occurrence) 데이터를 로딩한다.
# 이때 보통은 GBIF 에서 spocc 패키지를 이용해 받는다. (밑에 코멘트된 라인들을 이용해서)
# 하지만 이 실습에서는 미리 다운받은 데이터를 이용한다. 이 데이터는 패키지 기본데이터이기때문에
# 패키지를 로딩할때 이미 bvariegatus란 이름으로 불러져있다.
occs <- bvariegatus
# 발생 데이터가 어떻게 생겼는지 첫 몇줄을 보자. 발생지역의 좌표들이다.
head(occs)

# 같은 좌표를 가진 발생점(occurrence point) 들은 동일한 데이터가 복수되있는 것일 수 있기 때문에
# 없에주는게 좋은 습관이다.
occs <- occs[!duplicated(occs),]

# WorldClim의 기후 데이터를 사용하여 종의 기후적 적합성을 모델링해보자.
# WorldClim은 다양한 해상도에서 사용 가능한 여러 변수를 제공한다;
# 여기에서는 간단히 하기 위해 WorldClim 1.4의 dismo 패키지에 포함된 10arcmin 해상도
# (적도에서 약 20 km)의 9개 생물기후 변수를 사용해보자.
# 이 기후 데이터는 1950년부터 2000년까지의 50년 평균을 기반으로 하고 있다.

# 먼저 dismo 폴더에서 예측 변수 래스터 파일을 찾는다.
envs.files <- list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                         pattern='grd', full.names=TRUE)


# RasterStack에 래스터 파일을 읽는다.
# 이 변수들은 8개의 생물기후 변수와 하나의 범주형 변수 (categorical variable) "biome"을 나타낸다.
# 생물기후 변수의 설명은 다음에서 찾을 수 있다:
# https://www.worldclim.org/data/bioclim.html
envs <- raster::stack(envs.files)

# biome 래스터에는 다른 래스터에서 값이 있는 셀에 대해 NA가 몇 개 있다.
# 모든 래스터에 대해 이러한 셀의 값을 NA로 변경하기 위해 모든 래스터를 biome에 마스킹하자.
# ENMeval은 이를 자동으로 수행할 것이지만, 나중에 경고 메시지를 피하기 위해 여기에서 수행한다.
# RasterBrick에서 RasterStack으로 다시 변경하는 이유는 RasterBricks가 factor 래스터를
# 처리하지 못하기 때문이다.
envs <- raster::mask(envs, envs[[9]]) %>% raster::stack()

# 범주형 변수를 factor 로 선언하는 것을 확인한다.
envs$biome <- raster::as.factor(envs$biome)

# 이제 예측 변수 래스터에서 격자 셀을 공유하는 occurrence들을 제거하자.
# Maxent는 이를 기본적으로 수행하지만, 다른 알고리즘의 경우 연구의 목적에 따라
# 이 작업을 수행할지 말지를 결정해야 할 것이다.
# 공간적 자기상관을 (spatial autocorrelation) 피하기 위해 발생 기록을 서로 일정한 거리만큼 떨어뜨리는 또
# 다른 방법은 공간적 얇게 하기(spatial thinning)를 사용하는 것이다 (Aiello-Lammens et al. 2015).
occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]

# Rasterstack의 첫 번째 래스터, 연평균 온도를 그려보자.
plot(envs[[1]], main="Mean annual temperature")

# 래스터 위에 모든 발생 점을 추가한다.
points(occs)


# 아마존 강의 동쪽에 몇몇 점이 있다.
# 이것이 모델에 포함시키고 싶지 않은 개체군이라고 가정하자.
# 위도와 경도로 발생 사례를 부분집합으로 만들어 이러한 점들을 분석에서 제거할 수 있다.
occs <- filter(occs, latitude > -20, longitude < -45)

# 부분집합으로 만든 발생 사례를 그려 올바르게 필터링했는지 확인한다.
points(occs, col = 'red')

# 우리 모델은 발생 지역의 환경과 배경 지역의 환경을 비교할 것이므로, 
# 배경 범위에서 무작위 점을 샘플링해야 한다. 이제 글로벌 예측 변수 래스터를
# 더 작은 지역으로 자르는 것으로 연구 범위를 지정하고 배경 데이터를 샘플링할
# 위치를 정의할 것이다. 이동 제약과 같은 한계로 인해 적합하지만 빈번하지 않은
# 지역을 포함하지 않기 위해, 우리는 발생 지역을 둘러싼 지역으로 배경 범위를
# 보수적으로 정의할 것이다 (VanDerWal et al. 2009, Merow et al. 2013). 이를 위해
# 모든 발생 지역을 포함하는 경계 상자를 버퍼링할 것이다. 배경 범위를 구분하는 
# 다른 방법들(예: minimum convex hulls)은 점을 포함하는 지리적 공간을 더 잘 특성화하기
# 때문에 더 보수적이다. 어쨌든, 이것은 연구를 할 때 신중히 고려해야 할 많은
# 것들 중 하나이다.





# 이제 sf (simple features)라는 다른 공간 R 패키지를 사용할 것이다.
# occs를 sf 객체로 만들자 -- 이 점들의 좌표 참조 시스템(crs)은 WGS84로,
# 지리적 crs(위도/경도``)이며 우리의 envs 래스터의 crs와 동일하므로,
# 이를 RasterStack의 crs로 지정한다.
occs.sf <- sf::st_as_sf(occs, coords = c("longitude","latitude"), crs = raster::crs(envs))


# 이제, 점 데이터를 동일 면적 투영으로 변환한다. 
# 이는 우리의 도(degree)를 미터로 변환하는 것으로, 
# 버퍼링(다음 단계)에 이상적이# 우리는Eckert IV 투영을 사용한다.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

# 모든 발생 사례를 500 km로 버퍼링하고, 다각형을 합친다(시각화를 위해).
# 그리고 래스터 패키지가 사용할 수 있는 형태로 다시 변환한다. 마지막으로,
# 버퍼를 WGS84(위도/경도)로 다시 투영한다.
# 여기서는 카리브해 제도를 샘플링하지 않기 위해 500 km를 선택한다.
occs.buf <- sf::st_buffer(occs.sf, dist = 500000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(occs)

# sf 객체를 플롯에 추가하려면, add = TRUE를 사용한다.
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# 연구 범위와 일치하도록 환경 래스터를 자른다.
envs.bg <- raster::crop(envs, occs.buf)
# 다음으로, 버퍼의 형태에 맞게 래스터를 마스킹한다.
envs.bg <- raster::mask(envs.bg, occs.buf)
plot(envs.bg[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)


# 다음 단계에서는 배경에서 10,000개의 무작위 점을 샘플링할 것이다.
# (배경 점의 수도 연구에 따라 고려해야 할 사항임을 주의하라).

# 하나의 배경 범위 래스터에서 10,000개의 배경 점을 무작위로 샘플링한다(교체 없이 셀 당 하나만). 
# 주의: 래스터에 <10,000 픽셀이 있으므로 경고가 표시되고 모든 픽셀이 배경으로 사용될 것이다.
# 우리는 일부 격자 셀이 NA로되어 있기 때문에 바이옴 변수에서 샘플링할 것이며, NA를 
# 가진 배경 점을 얻는 것을 피하려고 한다. 스택의 한 래스터가 다른 래스터에 데이터가 있는
# 곳에서 NA를 가지고 있다면, ENMeval은 내부적으로 이러한 셀을 NA로 변환한다.
bg <- dismo::randomPoints(envs.bg[[9]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

# 배경점이 상당히 광범위하게 지역을 커버하고 있다는 것에 주목하라(모든 셀).
plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)

# 2. 모델 평가를 위한 발생 데이터 분할 (Partitioning Occurrences for Evaluation).

# 데이터 분할은 사용자가 partitions 인수에 대해 입력한 내용을 기반으로 
# 모델 튜닝과 (Model tuning) 모델 피팅 (model fitting)을 실행하는 ENMeval 패키지에서 가장 
# 주요 함수인 ENMevaluate()에 의해 내부적으로 수행되지만, 이것은 또한 분할목적을 위한
# 함수를 사용하여 외부적으로도 수행될 수 있다. 이 섹션에서는 이러한 다양한 함수를 
# 설명하고 보여준다.
# 또한 분할과 배경이나 연구 범위에 대한 환경적 유사성을 보여주는 정보
# 제공적인 그림을 만드는 방법을 보여준다.
# ENMeval 에서는 7가지 분할방법을 유치 한다.
# 1.Spatial Block
# 2.Spatial Checkerboard
# 3.Spatial Hierarchical Checkerboard
# 4.Jackknife (leave-one-out)
# 5.Random k-fold
# 6.Fully Withheld Testing Data
# 7.User
# 이 실습에선 1,2,5 번을 다뤄볼 것이다.

# a) Block
# 첫째로, 'block' 방법은 발생 지역을 위도와 경도 선에 따라 분할하여
# 가능한 한 동일한 수의 네 개의 공간 그룹으로 나눈다.
# 발생과 배경 지역 모두 이러한 선에 대한 위치를 기반으로 
# 네 개의 구간 각각에 할당된다 - 첫 번째 방향은 점들을 두 그룹으로 나누고,
# 두 번째 방향은 이것을 더욱 두 그룹으로 나눈다,
# 결과적으로 네 개의 그룹이 생성된다.
# 결과 객체는 각 발생과 배경 지점에 대한 구간 지정을 제공하는 두 개의 벡터의
# 목록이다. 
# ENMeval 2.0에서 사용자는 orientation 인수로 블로킹의
# 다른 방향을 추가로 지정할 수 있다.
block <- get.block(occs, bg, orientation = "lat_lon")


# 각 분할에 동일한 개수의 발생 데이터가 있는지 확인한다.
table(block$occs.grp)


# 우리는 우리의 예측 변수 래스터 중 하나에 분할을 플롯하여
# 그것들이 공간에서 어디에 떨어지는지 시각화할 수 있다.
# ENMeval 2.0의 그래프 함수는 ggplot2 (Wickham 2016)를 사용하며, 많은 온라인 자원이 있는
# R의 인기 있는 그래핑 패키지이다.
# 더하기 방식으로 다른 ggplot 함수를 ggplot에 추가할 수 있어,
# 이러한 쉽게 커스터마이징 할 수 있다.
evalplot.grps(pts = occs, pts.grp = block$occs.grp, envs = envs.bg) + 
  ggplot2::ggtitle("Spatial block partitions: occurrences")

# 배경을 그래프하면 배경 범위가 네 개의 구간에 걸쳐 점의 균일성을 극대화하는
# 방식으로 분할되어 있음을 보여준다. 면적의 균일성을 극대화하기 위한 것이 아니다.
evalplot.grps(pts = bg, pts.grp = block$bg.grp, envs = envs.bg) + 
  ggplot2::ggtitle("Spatial block partitions: background")

# 각 분할이 다른 모든 분할과 얼마나 다른 환경과 연관이 있는지 궁금하다면,
# evalplot 함수를 사용하여 각 분할을 참조로 하는 MESS 예측의 히스토그램이나 래스터를 플롯할 수 있다.
# 먼저 발생 지역과 배경 지역에서 예측 변수 값을 추출해야 한다.
# MESS 는 Multivariate Environmental Similarity surface 로 여러
# 변수들을 모두 고려해 셀마다 환경차이가 얼마나 나는지 알아볼 수 있게 해준다
# (Elith et al. 2010). 
occs.z <- cbind(occs, raster::extract(envs, occs))
bg.z <- cbind(bg, raster::extract(envs, bg))
evalplot.envSim.hist(sim.type = "mess", ref.data = "occs", occs.z = occs.z, bg.z = bg.z, 
                     occs.grp = block$occs.grp, bg.grp = block$bg.grp, categoricals = "biome")


# b) Checkerboard

# 다음 분할 방법은 발생 지역을 분할하기 위한 '체커보드' 접근법의 변형이다
# (Radosavljevic & Anderson 2014).
# 이 방법은 연구 범위에 체커보드 격자를 생성하고 체커보드에서 어디에 떨어지는지에 
# 따라 지역을 그룹으로 나눈다. 블록 방법과는 달리 두 체커보드 방법은 모두 지리적
# 공간을 동일하게 세분화하지만 각 구간에 발생 지역의 균형있는 수를 보장하지 않는다.
# 이 방법들에 대해 사용자는 기본 체커보드 패턴을 기반으로 하는 래스터 레이어를
# 제공해야 한다. 여기서는 단순히 예측 변수 RasterStack을 사용한다.
# 또한 사용자는 aggregation.factor를 정의해야 한다. 이 값은 기본 체커보드
# 패턴을 만들 때 집계할 그리드 셀의 수를 지정한다.

# 두가지 체커보드 분할 방식이 ENMeval에서 제공되는데 여기선 체커보드1만 보겠다.
# 방법은 간단한 체커보드 패턴을 사용하여 포인트를 k = 2 공간 그룹으로 분할한다.

cb1 <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor=5)
evalplot.grps(pts = occs, pts.grp = cb1$occs.grp, envs = envs.bg)

# 배경 포인트를 그래프하면 체커보드 패턴이 명확하게 보인다.
evalplot.grps(pts = bg, pts.grp = cb1$bg.grp, envs = envs.bg)

# MESS 맵에서 볼 수 있듯이 이 방법은 분할 간에 유사한 환경 표현을 결과로 가져온다.
evalplot.envSim.hist(sim.type = "mess", ref.data = "occs", occs.z = occs.z, bg.z = bg.z, 
                     occs.grp = cb1$occs.grp, bg.grp = cb1$bg.grp, categoricals = "biome")

# c) Random k-fold
# '랜덤 k-폴드' 방법은 사용자가 지정한 수의 (k) 빈에 발생 지역을 무작위로 분할한다
# (Hastie 등, 2009). 이 방법은 현재 버전의 Maxent 소프트웨어 GUI에서 사용 가능한 
# '크로스-밸리데이트' 분할법과 동일하다. 특히 더 큰 발생 데이터셋에서 이 분할
# 방법은 무작위로 일부 공간적 군집을 초래할 수 있으며, 이로 인해 공간 자기상관을 
# 처리하기 위해 공간 분할 방법이 선호된다(Roberts 등, 2017). 아래에서는 데이터를 
# 다섯 개의 무작위 그룹으로 분할한다.

rand <- get.randomkfold(occs, bg, k = 5)
evalplot.grps(pts = occs, pts.grp = rand$occs.grp, envs = envs.bg)

# 분할이 무작위이므로 분할마다 사이에 큰 환경적 차이는 없다.
evalplot.envSim.hist(sim.type = "mess", ref.data = "occs", occs.z = occs.z, 
                     bg.z = bg.z, occs.grp = rand$occs.grp, bg.grp = rand$bg.grp, 
                     categoricals = "biome")


# 3. ENMeval을 통한 Maxent 실행
# 여기에서는 Maxent 모델(maxent.jar 또는 maxnet) 실행과 튜닝 절차를 설명하겠습니다.
# Maxent 모델을 사용하기 위해 ENMevaluate를 호출할 때 조정해야
# 하는 두 가지 주요 매개변수가 있습니다.
# (1) 정규화 승수(regularization multiplier) 값의 범위와 
# (2) 특징 클래스의 (feature class) 조합입니다.
# 정규화 승수(RM)는 모델에 변수 또는 그 변형을 포함시키는 것에 대한 패널티를 결정합니다.
# 높은 RM 값은 모델의 복잡성에 더 강한 패널티를 부과하여 더 간단한(평탄한)
# 모델 예측을 만듭니다. 특징 클래스는 마진 반응 곡선의 (marginal response curves) 
# 잠재적인 모양을 결정합니다.
# 선형 특징 클래스만을 포함할 수 있는 모델은 모든 가능한 특징 클래스를 포함할 수 
# 있는 모델보다 더 간단할 가능성이 높습니다.

# ENMevaluate()는 RM 값과 특징 클래스 조합 각각에 대해 별도의 모델을 구축한다.
# 예를 들어, 다음 호출은 2개의 모델을 구축하고 평가한다. 하나는 RM=1을 사용하고
# 다른 하나는 RM=2를 사용하며, 둘 다 선형 특징만을 허용한다.
e.mx.l <- ENMevaluate(occs = occs, envs = envs, bg = bg, 
                      algorithm = 'maxnet', partitions = 'block', 
                      tune.args = list(fc = "L", rm = 1:2))

# occs 파라미터에는 발생 데이터, envs 에는 환경 레스터 데이터,bg 에는 배경 데이터를 입력해준다.
# 여기서 algorithm 이란 파라미터에는 'maxnet' 나 'maxent.jar' string 값을 넣어
# Maxent 를 사용할 수 있다.
# Maxent.jar 를 사용하기 위해선 다음 라인을 통해 rmaxent패키지를 받아야한다. (커맨트 처리되어있음)
# devtools::install_github("johnbaums/rmaxent")

# partitions parameter 에는 어떤 데이터 분할 방법을 택할지 정한다. 
# 공간적이지 않은 분할 방법으론 "randomkfold", "jackknife" 가있고
# 공간적 분할 방법으론 "block", "checkerboard1", "checkerboard2", 가있고
# "testing" 은 fully withheld data, "user" 는 유저 커스터마이징, "none" 는 분할없음이다. 

# 더 다양한 특성 클래스와 정규화 멀티플라이어를 사용할 수 있는 
# 더 넓은 범위의 모델을 비교하고 싶을 수도 있다. 밑에 줄은 4 가지의 특징 클레스와 5개의 RM값을 
# 조합으로 모델 튜닝을 한다. 
# 밑 커멘트 된 줄은 실행하는데 너무 오래 걸리니, 미리 실행해 놓은 모델 객체를 불러온다.
#e.mx <- ENMevaluate(occs = occs, envs = envs, bg = bg, 
#                    algorithm = 'maxnet', partitions = 'block', 
#                    tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5))
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
e.mx <- read_rds('./data/emx.rds')

# 모델 튜닝 결과를 보자.
# 다음은교차 검증을 위한 테스트 데이터에 대한 요약 통계를 포함한 결과 테이블의 일부다.
eval.results(e.mx) %>% head()

# 각 열은 그 열에 있는 매개변수 조합 (tune.args) 을 써서 MaxEnt알고리즘을 실행해 얻은 모델에 대한
# 평가 통계를 포함한 결과 테이블이다.
# AUC, AUCdiff, CBI, OR_10, OR_mtp, AICc 등 다양한 평가 통계 척도들의 있는
# 것을 볼 수 있다.


# ggplot의 facetting 함수를 사용하여 여러 통계를 그래프할 수 있다.
evalplot.stats(e = e.mx, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm")


# 가끔 오류 막대가 그래프를 보기 어렵게 만들 수 있으므로, 이를 끄는 것을 시도할 수 있다.
evalplot.stats(e = e.mx, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               error.bars = FALSE)



# 4. 모델 선택 (Model selection)

# 이제 실행한 다른 매개변수 조합으로 만든 모델 중에서 하나 이상의 최적의 모델을 선택해야 한다.
# 여러 매개변수 조합으로 모델을 만들어보고 그중에서 최적의 모델을 고름으로서 주어진 데이터에 대한
# 최적의 매개변수 조합을 찾는 것이 모델 튜닝인 것이다.
# 이 예제에서는 두가지 모델 선택 방법을 보겠다.
# 1) AICc를 사용하여 교차 검증 결과를 고려하지 않고 모델을 선택하는 방법과 
# 2) 평균 Omission Rate 이 가장 낮고, 평균 AUC가 가장 높은 모델을 선택하는 
# 순차적 방법을 보여줄 것이다.
# 검증 AUC가 절대적인 성능 측정 지표로서는 부적절하다고 지적된바 있지만,
# 동일한 데이터로 구성된 모델 간의 상대적인 비교에는 유효하다. 
# 또 ENMeval 2.0도 R 패키지 ecospat을 사용하여 
# Continuous Boyce Index (CBI), 또는 줄여서 Boyce Index,
# 통계지표를 계산하기 때문에 이 지표도 최적의 모델을 선택하는 데 사용할 수 있다.

# 결과 테이블을 개체로 지정한다.
res <- eval.results(e.mx)

# 방법 1)
# delta AICc 점수가 0인 모델을 선택하거나 AICc 점수가 가장 낮은 모델을 선택한다.
# 실제로는 AICc 점수의 차이가 2 미만인 모델은 일반적으로 통계적으로 동등하다고 간주된다.
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

# 방법 2)
# 이 줄은 위에서 설명한 순차적 기준을 실행한다.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

# 이제 순차적 기준(방법2)을 기반으로 최적의 모델 설정을 선택하고 살펴본다.
# 최적의 모델의 tune.args를 사용하여 ENMevaluation 객체에서 단일 모델을 선택한다.
mod.seq <- eval.models(e.mx)[[opt.seq$tune.args]]
# 모델에서 0이 아닌 계수를 본다.
mod.seq$betas

# 모델에서 0이 아닌 계수를 가진 예측 변수들에 대한 Marginal Response Curves다.
# y축을 cloglog 아웃풋으로 정의하는데, 이는 0과 1 사이에 제한된 출현 확률의 
# 근사값이다.
plot(mod.seq, type = "cloglog")
dev.off()



# 우리는 위에서 모델 객체에 대해 한 것과 동일한 방식으로 최적 모델의 예측을 
# 선택할 수 있다.
pred.seq <- eval.predictions(e.mx)[[opt.seq$tune.args]]
plot(pred.seq)

# 우리는 훈련 데이터가 어디에 위치하는지 시각화하기 위해 발생 지점 
# 위에 배경 지점을 구간별로 표시할 수도 있다.
points(eval.occs(e.mx), pch = 21, bg = eval.occs.grp(e.mx))
points(eval.bg(e.mx), pch = 3, col = eval.bg.grp(e.mx), cex = 0.5)
