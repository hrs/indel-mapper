# Indel Mapper

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

## Set up

First, make sure the correct version of `virtualenv` is installed:

```shell
$ pip3 install virtualenv
```

Next, `cd` into your project and set up the `virtualenv` directory:

```shell
$ virtualenv --python=python3 indel-mapper
```

Activate `virtualenv`. This adds the `indel-mapper/bin` directory to the start
of your `$PATH`.

```shell
$ source indel-mapper/bin/activate
```

Install the required libraries:

```shell
$ pip3 install -r requirements.txt
```

Install the Redis server. The instructions for this will vary depending on your
operating system. For Debian, Ubuntu, and derivatives,

```shell
$ sudo apt install redis-server
```

## Starting the Web application

Run the application:

```shell
$ python3 app.py
```

Start Redis. How you do this depends on your operating system.

**Debian, Ubuntu, or derivatives:**

```shell
$ sudo systemctl start redis-server
```

Run Celery

```shell
$ indel-mapper/bin/celery worker -A app.celery --loglevel=info
```

## Running the test suite

Run the tests:

```shell
$ python3 -m unittest
```

## Deploying to Heroku

```shell
$ git push heroku master
```

## Running locally as a command line application

Example:

```shell
$ python3 run.py -a ~/Documents/bowtie2_results.sam -r ~/Documents/references.csv -o ~/Documents/results.csv
```
There are three required arguments:

* `-a` or `--alignment` Alignment SAM file
* `-r` or `--reference` Reference CSV file
* `-o` or `--output` Output file, in CSV

## License

Indel Mapper is licensed under Version 3 of the GNU General Public License.
