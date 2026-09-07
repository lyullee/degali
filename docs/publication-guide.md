# DEGALI 공개 배포 가이드: GitHub, DOI, PyPI

기준일: 2026-09-07. 저장소는 `lyullee/degali`, 패키지명은 `degali`다.

## 1. 권장 공개 순서

```text
GitHub 저장소 생성과 alpha 소스 공개
  -> Zenodo 초안에서 DOI 선예약
  -> DOI·정식 버전·배포일을 메타데이터에 반영
  -> v0.1.0 커밋·태그·GitHub Release
  -> PyPI Trusted Publishing 자동 배포
  -> 같은 산출물을 Zenodo 초안에 올리고 DOI 발행
  -> 공개 페이지와 해시 최종 확인
```

이 순서는 첫 배포물 자체에 DOI를 넣을 수 있고, 장기 API 토큰을 저장하지
않는다. Zenodo의 GitHub 자동보관을 사용해 먼저 릴리스하면 DOI는 릴리스 뒤에
생기므로 첫 패키지의 `CITATION.cff`에 그 DOI가 들어가지 않는다. 첫 버전은
수동 DOI 선예약 방식을 권장하며, 같은 릴리스를 GitHub 자동보관으로 다시
등록해 DOI를 중복 발행하지 않는다.

## 2. GitHub 저장소

1. <https://github.com/new>에서 저장소를 만든다.
2. Owner는 `lyullee`, Repository name은 `degali`, Visibility는 `Public`으로 한다.
3. README, `.gitignore`, license 자동 생성을 모두 끈다. 로컬에 이미 있다.
4. 생성 후 로컬에서 다음을 실행한다.

```bash
git commit -m "Initial public alpha snapshot of DEGALI"
git remote add origin https://github.com/lyullee/degali.git
git push -u origin main
```

5. GitHub의 Actions 탭에서 공개 데이터가 필요 없는 `tests`가 통과하는지 확인한다. 로컬 Fortran
   오라클을 포함하므로 로컬 Windows 시험보다 이 결과가 공개 검증의 기준이다.
6. Settings > Environments에서 `pypi` 환경을 만들고, 가능하면 본인을
   required reviewer로 지정한다. 공개 릴리스가 곧바로 PyPI에 올라가는 것을
   막는 마지막 승인 장치다.

## 3. Zenodo DOI 선예약

1. <https://zenodo.org>에 GitHub 또는 ORCID로 로그인한다.
2. Profile > Linked accounts에서 GitHub를 연결한다.
3. 새 upload를 만들고 Upload type을 `Software`로 선택한다.
4. 제목은 `DEGALI: Dense Gas Dispersion for Liquid Hydrogen`로 입력한다.
5. Creator는 `Lee, Ugwiyeon`, affiliation은
   `Korea Gas Safety Corporation`으로 입력한다.
6. License는 `MIT`, access는 `Open`, version은 `0.1.0`으로 둔다.
7. DOI 항목에서 기존 DOI가 없다고 선택하고 `Get a DOI now!`로 DOI를
   선예약한다. 이 단계에서는 아직 Publish를 누르지 않는다.
8. 발급된 DOI를 정확히 복사해 다음 위치에 넣는다.

```text
CITATION.cff: doi, version, date-released
pyproject.toml: project.urls의 DOI 또는 Archive 항목
README.md: Zenodo DOI badge와 인용 링크
```

9. 정식 배포물 생성 후 아래 두 파일을 같은 Zenodo 초안에 업로드한다.

```text
dist/degali-0.1.0.tar.gz
dist/degali-0.1.0-py3-none-any.whl
```

필요하면 GitHub가 생성한 source archive도 함께 올릴 수 있으나 같은 내용을
불필요하게 중복하지 않는다. 파일과 메타데이터를 확인한 뒤에만 Zenodo
`Publish`를 누른다. Zenodo 공개 후에는 그 버전의 파일을 교체할 수 없다.

## 4. PyPI Trusted Publishing

PyPI API token을 로컬이나 GitHub secret에 저장하지 않는다. 저장소의
`.github/workflows/release.yml`은 GitHub OIDC와 PyPI Trusted Publishing을
사용하며, 배포 증명(attestation)도 기본 생성한다.

1. <https://pypi.org/account/register/>에서 계정을 만들거나 로그인하고
   이메일를 활성화한다.
2. Account settings > Publishing > Add a new pending publisher로 이동한다.
3. 다음 값을 정확히 입력한다.

| PyPI 필드 | 값 |
|---|---|
| PyPI project name | `degali` |
| Owner | `lyullee` |
| Repository name | `degali` |
| Workflow filename | `release.yml` |
| Environment name | `pypi` |

4. Pending publisher는 이름을 예약하지 않는다. 설정 후 첫 배포를 너무 오래
   미루지 않는다.
5. 첫 업로드가 성공하면 pending publisher가 일반 trusted publisher로
   자동 전환되는지 PyPI의 Publishing 설정에서 확인한다.

## 5. 첫 정식 버전 준비

DOI를 선예약한 뒤 개발 버전을 정식 버전으로 바꾼다.

```text
pyproject.toml: version = "0.1.0"
CITATION.cff: version = "0.1.0"
CITATION.cff: date-released = 2026-09-07
CITATION.cff: doi = "선예약 DOI"
```

그다음 전체 공개 검사와 핵심 시험을 실행한다.

```bash
python tools/check_publication.py
python -m pytest -m "not slow" -q
python -m build
python -m twine check dist/*
```

정식 커밋과 태그는 버전과 DOI가 들어간 동일 상태에서 만든다.

```bash
git add .
git commit -m "Release DEGALI 0.1.0"
git tag -a v0.1.0 -m "DEGALI 0.1.0"
git push origin main
git push origin v0.1.0
```

## 6. GitHub Release와 자동 PyPI 배포

1. GitHub > Releases > Draft a new release를 연다.
2. 기존 태그 `v0.1.0`을 선택한다. 새로 다른 태그를 만들지 않는다.
3. 제목은 `DEGALI 0.1.0`, 설명에는 alpha 연구용 경고와 CHANGELOG를 넣는다.
4. Publish release를 누르면 `release.yml`이 실행된다.
5. workflow는 태그와 패키지 버전이 정확히 같은지, `.dev0`가 아닌지,
   공개 검사와 `twine check`가 통과하는지 확인한 뒤 PyPI에 올린다.
6. `pypi` 환경에 reviewer를 지정했다면 Actions 화면에서 승인한다.
7. <https://pypi.org/project/degali/>에서 파일, 저자, 프로젝트 링크와
   provenance/attestation을 확인한다.

## 7. 마지막 일치 검사

다음 네 위치가 모두 `0.1.0`, 동일 DOI와 동일 저장소를 가리켜야 한다.

- GitHub Release `v0.1.0`
- PyPI `degali 0.1.0`
- Zenodo software record `0.1.0`
- 저장소의 `CITATION.cff`

로컬 `dist` 파일과 공개 파일의 SHA-256도 비교한다. 문제가 생기면 기존 PyPI
버전과 Zenodo 파일을 덮어쓰지 말고 코드를 수정해 `0.1.1`을 발행한다.

## 8. 현재 중단점

- GitHub CLI는 현재 PC에 설치돼 있지 않다. 웹에서 저장소를 만들거나 CLI를
  설치한 뒤 인증해야 한다.
- `lyullee` 계정은 Ugwiyeon Lee의 소유로 확인됐다.
- Zenodo 계정 로그인과 DOI 선예약은 사용자 본인 확인이 필요하다.
- PyPI 계정 로그인, 2FA 및 pending trusted publisher 등록은 사용자 본인
  확인이 필요하다.
- DOI가 아직 없으므로 버전은 의도적으로 `0.1.0.dev0`로 유지한다.
