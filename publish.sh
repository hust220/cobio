#python setup.py sdist bdist_wheel
python setup.py bdist_wheel

for whl in dist/*.whl; do
  auditwheel repair "$whl" -w dist
  rm "$whl"
done

#twine upload dist/*
