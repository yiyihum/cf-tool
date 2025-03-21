package cmd

import (
	"io/ioutil"

	"github.com/yiyihum/cf-tool/client"
	"github.com/yiyihum/cf-tool/config"
)

// Submit command
func Submit() (err error) {
	cln := client.Instance
	cfg := config.Instance
	info := Args.Info
	submitOnly := Args.SubmitOnly
	filename, index, err := getOneCode(Args.File, cfg.Template)
	if err != nil {
		return
	}

	bytes, err := ioutil.ReadFile(filename)
	if err != nil {
		return
	}
	source := string(bytes)

	lang := cfg.Template[index].Lang
	if err = cln.Submit(info, lang, source, submitOnly); err != nil {
		if err = loginAgain(cln, err); err == nil {
			err = cln.Submit(info, lang, source, submitOnly)
		}
	}
	return
}
